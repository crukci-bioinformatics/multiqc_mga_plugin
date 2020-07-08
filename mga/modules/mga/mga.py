#!/usr/bin/env python

""" MultiQC Multi Genome Alignment module """

import colorsys
import logging
import operator
import os

from collections import OrderedDict
from lxml import objectify, etree
from types import SimpleNamespace

from multiqc import config
from multiqc.plots import bargraph, table
from multiqc.modules.base_module import BaseMultiqcModule

from .colour import Colour

log = config.logger

# See https://stackoverflow.com/a/16279578
bar_colours = SimpleNamespace(**{
    'red': Colour.fromBytes(255, 0, 0),
    'orange': Colour.fromBytes(255, 200, 0),
    'green': Colour.fromBytes(0, 255, 0),
    'grey': Colour.fromBytes(128, 128, 128),
    'adapter': Colour.fromBytes(255, 102, 255)
})

max_alpha = 1.0
min_alpha = 0.25
max_error = 0.0025
min_error = 0.01

assigned_fraction_threshold = 0.01
aligned_fraction_threshold = 0.01
error_rate_threshold = 0.0125

# Based on https://github.com/MultiQC/example-plugin

class MultiqcModule(BaseMultiqcModule):
    
    def __init__(self):

        # Halt execution if we've disabled the plugin
        if config.kwargs.get('disable_plugin', True):
            return None

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name = 'Multi Genome Alignment',
            target = "MGA",
            anchor = 'mga',
            href = "https://github.com/crukci-bioinformatics/MGA",
            info = "is used to align a sample of reads to multiple genomes for contamination screening."
        )
        
        self.mga_data = dict()
        for mgafile in self.find_log_files('mga', filecontents=False, filehandles=True):
            self.read_mga_xml_file(mgafile)

        if len(self.mga_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.mga_data)))
        
        sum_sequences = etree.XPath("sum(//MultiGenomeAlignmentSummaries/MultiGenomeAlignmentSummary/SequenceCount)")
        count_references = etree.XPath("count(//MultiGenomeAlignmentSummaries/ReferenceGenomes/ReferenceGenome)")
        
        for run_id, dataset in self.mga_data.items():
            dataset_props = self._read_properties(dataset)
            
            total_sequence_count = sum_sequences(dataset)
            yield_multiplier = 2 if dataset_props.get('endtype') == 'Paired End' else 1
            cycles = int(dataset_props['cycles'])
            total_yield = yield_multiplier * cycles * total_sequence_count / 1000000000.0
            
            bar_data, bar_categories, plot_config = self.get_bar_data(dataset)

            bargraph_plot_html = bargraph.plot(bar_data, bar_categories, plot_config)
            
            trim_start = int(float(dataset.findtext("TrimStart")))
            trim_length = int(float(dataset.findtext("TrimLength")))
            number_of_genomes = int(count_references(dataset))

            self.add_section(
                name = 'Sequencing Results',
                description = '''
                    <table>
                        <tr><td>Flow Cell ID:</td><td>{}</td></tr>
                        <tr><td>Run name:</td><td>{}</td></tr>
                        <tr><td>Cycles:</td><td>{}</td></tr>
                        <tr><td>End type:</td><td>{}</td></tr>
                        <tr><td>Yield (Gbases):</td><td>{:.2f}</td></tr>
                        <tr><td>Total sequences:</td><td>{:,.0f}</td></tr>
                    </table>
                    <br/>
                '''.format(dataset_props.get('flowcellid'), dataset_props.get('runname'), dataset_props.get('cycles'),
                           dataset_props.get('endtype'), total_yield, total_sequence_count),
                anchor = 'mga-plots',
                helptext = '''
                    Sequences were sampled, trimmed to {} bases starting from position {}, and mapped to {} reference genomes
                    (see list below) using Bowtie. Sequences containing adapters were found by ungapped alignment of the full
                    length sequence to a set of known adapter and primer sequences using Exonerate. Further details on the
                    alignment results and the assignment of reads to genomes are given below.
                '''.format(trim_length, trim_start, number_of_genomes),
                plot = bargraph_plot_html
            )

            for summary in dataset.findall("MultiGenomeAlignmentSummary"):
                dataset_id = summary.findtext("DatasetId")
                self.add_section(
                    name = 'Dataset "{}" Statistics'.format(dataset_id),
                    description = "Statistics for dataset {}".format(dataset_id),
                    anchor = "mga-stats-{}".format(dataset_id),
                    plot = self.get_table_data(summary)
                )
                
                self.add_section(
                    name = 'Dataset "{}" Samples'.format(dataset_id),
                    description = "Sample Details for dataset {}".format(dataset_id),
                    anchor = "mga-samples-{}".format(dataset_id),
                    plot = self.get_sample_table_data(summary)
                )


    def read_mga_xml_file(self, mgafile):
        content = objectify.parse(mgafile["f"])
        if content is None or content.getroot() is None:
            log.warn("Failed to parse {}".format(mgafile['fn']))
        elif content.getroot().tag != 'MultiGenomeAlignmentSummaries':
            log.debug("{}/{} is not an MGA summary file.".format(mgafile['root'], mgafile['fn']))
        else:
            log.debug("Found file {}/{}".format(mgafile["root"], mgafile["fn"]))
            run_id = content.findtext("RunID")
            if run_id not in self.mga_data:
                self.mga_data[run_id] = content


    def get_bar_data(self, mga_data):
        
        run_id = mga_data.findtext("RunID")
        bar_data = OrderedDict()
        categories = OrderedDict()
        max_sequenced_count = 0
        
        for summary in mga_data.findall("/MultiGenomeAlignmentSummary"):
            dataset_id = summary.findtext("DatasetId")
            sequence_count = int(summary.findtext("SequenceCount"))
            sampled_count = int(summary.findtext("SampledCount"))
            adapter_count = int(summary.findtext("AdapterCount"))
            unmapped_count = int(summary.findtext("UnmappedCount"))
            
            sampled_to_sequenced = sequence_count / sampled_count
            
            max_sequenced_count = max(max_sequenced_count, sequence_count)

            species, controls = self._get_species_and_controls(summary)
            
            log.debug("Dataset {} Species: {}".format(dataset_id, species))
            log.debug("Dataset {} Controls: {}".format(dataset_id, controls))
            
            dataset_bar_data = dict()
            dataset_categories = dict()
            
            for alignment_summary in summary.findall("AlignmentSummaries/AlignmentSummary"):
                if self._accept_genome(species, summary, alignment_summary):
                    aligned_count = int(alignment_summary.findtext("AlignedCount"))
                    aligned_error = float(alignment_summary.findtext("ErrorRate"))
                    assigned_count = int(alignment_summary.findtext("AssignedCount"))
                    assigned_error = float(alignment_summary.findtext("AssignedErrorRate"))
                    
                    reference_genome_id = alignment_summary.find("ReferenceGenome").attrib['id']
                    reference_genome_name = alignment_summary.find("ReferenceGenome").attrib['name']
    
                    category_id = "{}.{}".format(dataset_id, reference_genome_id)
                    
                    colour = bar_colours.red
                    if reference_genome_name in controls:
                        colour = bar_colours.orange
                    elif reference_genome_name in species:
                        colour = bar_colours.green
                    elif len(species) == 0 or 'other' in map(lambda s: s.lower(), species):
                        colour = bar_colours.grey
                    
                    alpha = max_alpha - (max_alpha - min_alpha) * (assigned_error - min_error) / (max_error - min_error)
                    alpha = 1.0 - max(min_alpha, min(max_alpha, alpha))
                    
                    if assigned_count >= 100:
                        log.debug("{}\t{}\t{}\t{}".format(reference_genome_id, assigned_count, aligned_error * 100.0, alpha))
    
                    dataset_bar_data[category_id] = assigned_count * sampled_to_sequenced
                
                    dataset_categories[category_id] = {
                        'name': reference_genome_name,
                        'color': colour.applyAlpha(alpha).toHtml()
                    }
                
            # Sort into decreasing order of assigned count.
            dataset_bar_data = OrderedDict(sorted(dataset_bar_data.items(), key = lambda x: -x[1]))
            
            # The order of categories matters for the order in the plot, so add them to the
            # overall categories dictionary in the order they are in for the bar data.
            for category_id in dataset_bar_data.keys():
                categories[category_id] = dataset_categories[category_id]
                    
            log.debug("Adapter count: {} / {}".format(adapter_count, sampled_count))
            
            dataset_bar_data['adapter'] = adapter_count * sampled_to_sequenced

            bar_data[dataset_id] = dataset_bar_data

        categories['adapter'] = {
            'name': 'Adapter',
            'color': bar_colours.adapter.toHtml()
        }
        
        log.debug("Bar data = {}".format(bar_data))
        log.debug("Categories = {}".format(categories))

        plot_config = {
            'id': "mga_plot_{}".format(run_id),
            'title': "Multi Genome Alignment: {}".format(run_id),
            'cpswitch_counts_label': 'Read Counts',
            'xlab': "Lanes",
            'ylab': "Reads",
            'ymin': 0,
            'ymax': max_sequenced_count,
            'use_legend': False,
            'tt_percentages': False
        }
        
        return bar_data, categories, plot_config

    def get_table_data(self, summary):
        headers = OrderedDict()
        headers['species'] = {
            'title': 'Species/Reference Genome',
            'description': 'Reference genome species',
            'scale': False
        }
        headers['aligned_count'] = {
            'title': 'Aligned',
            'description': 'Number of reads aligned',
            'min': 0,
            'format': '{:d}',
            'scale': False
        }
        headers['aligned_perc'] = {
            'title': 'Aligned %',
            'description': 'Percentage of reads aligned',
            'suffix': '%',
            'min': 0,
            'max': 100,
            'format': '{:,.1f}',
            'scale': False
        }
        headers['aligned_error'] = {
            'title': 'Error rate',
            'description': 'Aligned error rate',
            'suffix': '%',
            'min': 0,
            'max': 100,
            'format': '{:,.2f}',
            'scale': False
        }
        headers['unique_count'] = {
            'title': 'Unique',
            'description': 'Number of uniquely aligned reads',
            'min': 0,
            'format': '{:d}',
            'scale': False
        }
        headers['unique_error'] = {
            'title': 'Error rate',
            'description': 'Uniquely aligned error rate',
            'suffix': '%',
            'min': 0,
            'max': 100,
            'format': '{:,.2f}',
            'scale': False
        }
        headers['preferred_count'] = {
            'title': 'Best',
            'description': 'Number of reads for which the genome is the best alignment',
            'min': 0,
            'format': '{:d}',
            'scale': False
        }
        headers['preferred_error'] = {
            'title': 'Error rate',
            'description': 'Error rate for best aligned reads',
            'suffix': '%',
            'min': 0,
            'max': 100,
            'format': '{:,.2f}',
            'scale': False
        }
        headers['assigned_count'] = {
            'title': 'Assigned',
            'description': 'Number of reads assigned',
            'min': 0,
            'format': '{:d}',
            'scale': False
        }
        headers['assigned_perc'] = {
            'title': 'Assigned %',
            'description': 'Percentage of reads assigned',
            'suffix': '%',
            'min': 0,
            'max': 100,
            'format': '{:,.1f}',
            'scale': False
        }
        headers['assigned_error'] = {
            'title': 'Error rate',
            'description': 'Assigned error rate',
            'suffix': '%',
            'min': 0,
            'max': 100,
            'format': '{:,.2f}',
            'scale': False
        }

        dataset_id = summary.findtext('DatasetId')

        table_data = dict()
        
        sequence_count = int(summary.findtext('SequenceCount'))
        sampled_count = int(summary.findtext('SampledCount'))
        adapter_count = int(summary.findtext('AdapterCount'))
        unmapped_count = int(summary.findtext('UnmappedCount'))

        species, controls = self._get_species_and_controls(summary)
            
        number_of_others = 0
        other_assigned_count = 0
        
        for alignment_summary in summary.findall("AlignmentSummaries/AlignmentSummary"):
            if not self._accept_genome(species, summary, alignment_summary):
                assigned_count = int(alignment_summary.findtext("AssignedCount"))
                number_of_others = number_of_others + 1
                other_assigned_count = other_assigned_count + assigned_count
        
        for alignment_summary in summary.findall("AlignmentSummaries/AlignmentSummary"):
            if number_of_others < 2 or self._accept_genome(species, summary, alignment_summary):
                aligned_count = int(alignment_summary.findtext("AlignedCount"))
                aligned_error = float(alignment_summary.findtext("ErrorRate"))
                unique_count = int(alignment_summary.findtext("UniquelyAlignedCount"))
                unique_error = float(alignment_summary.findtext("UniquelyAlignedErrorRate"))
                preferred_count = int(alignment_summary.findtext("PreferentiallyAlignedCount"))
                preferred_error = float(alignment_summary.findtext("PreferentiallyAlignedErrorRate"))
                assigned_count = int(alignment_summary.findtext("AssignedCount"))
                assigned_error = float(alignment_summary.findtext("AssignedErrorRate"))
                
                reference_genome_id = alignment_summary.find("ReferenceGenome").attrib['id']
                reference_genome_name = alignment_summary.find("ReferenceGenome").attrib['name']
            
                aligned_fraction = float(aligned_count) / float(sampled_count)
                assigned_fraction = float(assigned_count) / float(sampled_count)
                
                table_data[reference_genome_id] = {
                    'species': reference_genome_name,
                    'aligned_count': aligned_count,
                    'aligned_perc': aligned_fraction * 100.0,
                    'aligned_error': aligned_error * 100.0,
                    'unique_count': unique_count,
                    'unique_error': unique_error * 100.0,
                    'preferred_count': preferred_count,
                    'preferred_error': preferred_error * 100.0,
                    'assigned_count': assigned_count,
                    'assigned_perc': assigned_fraction * 100.0,
                    'assigned_error': assigned_error * 100.0
                }
        
        # Sort into decreasing order of assigned count.
        table_data = OrderedDict(sorted(table_data.items(), key = lambda x: -x[1]['assigned_count']))
        
        if number_of_others >= 2:
            table_data['Other'] = {
                'species': "{} others".format(number_of_others),
                'aligned_count': other_assigned_count,
                'aligned_perc': other_assigned_count * 100.0 / sampled_count
            }
        
        table_data['Unmapped'] = {
            'species': '',
            'aligned_count': unmapped_count,
            'aligned_perc': unmapped_count * 100.0 / sampled_count
        }
        
        table_data['Adapter'] = {
            'species': '',
            'aligned_count': adapter_count,
            'aligned_perc': adapter_count * 100.0 / sampled_count
        }
        
        table_config = {
            'namespace': 'mga',
            'id': 'mga_stats_table.{}'.format(dataset_id),
            'table_title': dataset_id,
            'col1_header': 'Reference ID',
            'no_beeswarm': True,
            'sortRows': False
        }
        
        return table.plot(table_data, headers, table_config)
    
    
    def get_sample_table_data(self, summary):
        headers = OrderedDict()
        """
        headers['sample_id'] = {
            'title': 'Sample ID',
            'description': 'Sample identifier',
            'scale': False
        }
        headers['sample_name'] = {
            'title': 'Sample Name',
            'description': 'Sample name',
            'scale': False
        }
        """
        headers['group'] = {
            'title': 'Group',
            'description': 'Research group name',
            'scale': False
        }
        headers['owner'] = {
            'title': 'Owner',
            'description': 'Researcher name',
            'scale': False
        }
        headers['sequence_type'] = {
            'title': 'Sequence Type',
            'description': 'The type of material in the sample',
            'scale': False
        }
        headers['end_type'] = {
            'title': 'End Type',
            'description': 'The sequencing method: single read or paired end',
            'scale': False
        }
        headers['library_type'] = {
            'title': 'Experiment Type',
            'description': 'The type of library (technique) used to create the sequencing pool',
            'scale': False
        }
        headers['species'] = {
            'title': 'Species',
            'description': 'Reference genome species',
            'scale': False
        }
        headers['control'] = {
            'title': 'Control',
            'description': 'Whether the reference genome is a control reference or not',
            'scale': False
        }

        dataset_id = summary.findtext('DatasetId')

        table_data = OrderedDict()
        
        samples = self._get_sample_information(summary)

        for sample in samples:
            sample_id = sample['sampleid']
            sample_name = sample['samplename']
            row_id = "{} / {}".format(sample_id, sample_name)
            
            table_data[row_id] = {
                'group': sample.get('group'),
                'owner': sample.get('owner'),
                'sequence_type': sample.get('sequencetype'),
                'end_type': sample.get('endtype'),
                'library_type': sample.get('experimenttype'),
                'species': sample.get('species'),
                'control': 'Yes' if sample.get('control') else 'No'
            }
        
        table_config = {
            'namespace': 'mga',
            'id': 'mga_sample_table.{}'.format(dataset_id),
            'table_title': dataset_id,
            'col1_header': 'Sample ID and name',
            'no_beeswarm': True,
            'sortRows': False
        }
        
        return table.plot(table_data, headers, table_config)
    

    def _get_species_and_controls(self, summary):
        samples = self._get_sample_information(summary)
        
        species = set()
        controls = set()
        for sample in samples:
            if 'species' in sample:
                sp = sample['species']
                species.add(sp)
                if sample['control']:
                    controls.add(sp)
                    
        return species, controls

    
    def _accept_genome(self, species_set, summary, alignment_summary):
        reference_genome_name = alignment_summary.find("ReferenceGenome").attrib['name']

        if reference_genome_name in species_set:
            return True

        sampled_count = int(summary.findtext("SampledCount"))
        
        aligned_count = int(alignment_summary.findtext("AlignedCount"))
        aligned_error_rate = float(alignment_summary.findtext("ErrorRate"))
        assigned_count = int(alignment_summary.findtext("AssignedCount"))
        
        aligned_fraction = float(aligned_count) / float(sampled_count)
        assigned_fraction = float(assigned_count) / float(sampled_count)

        return assigned_fraction >= assigned_fraction_threshold or aligned_fraction >= aligned_fraction_threshold and aligned_error_rate < error_rate_threshold


    def _get_sample_information(self, summary):
        samples = []
        for sample in summary.findall("Samples/Sample"):
            sample_info = { 'control': 'No' }
            self._read_properties(sample, sample_info)
            sample_info['control'] = sample_info['control'].lower() in ['yes','true']
            samples.append(sample_info)
        return samples

    
    def _read_properties(self, element, props = dict()):
        for prop in element.findall("Properties/Property"):
            try:
                props[prop.attrib['name'].replace(' ', '').lower()] = prop.attrib['value']
            except KeyError:
                # Leave put any that have no value.
                pass
        return props


#!/usr/bin/env python

""" MultiQC Multi Genome Alignment module """

import colorsys
import logging
import operator
import os

from collections import OrderedDict
from lxml import objectify
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
min_alpha = 0.1
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
            name = 'Multi Genome Alignemnt',
            target = "mga",
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
        
        for run_id, dataset in self.mga_data.items():
            bar_data, bar_categories = self.get_bar_data(dataset)

            plot_config = {
                'id': "mga_plot_{}".format(run_id),
                'title': "Multi Genome Alignment: {}".format(run_id),
                'xlab': "Lanes",
                'ylab': "Reads",
                'use_legend': False
            }
            
            bargraph_plot_html = bargraph.plot(bar_data, bar_categories, plot_config)
            
            self.add_section(
                description = 'MGA plots',
                helptext = '''
                Copy from the existing report.
                ''',
                plot = bargraph_plot_html
            )

            for summary in dataset.findall("MultiGenomeAlignmentSummary"):
                self.add_section(
                    description = "Statistics for dataset {}".format(summary.findtext("DatasetId")),
                    plot = self.get_table_data(summary)
                )


    def read_mga_xml_file(self, mgafile):
        content = objectify.parse(mgafile["f"])
        if content is None or content.getroot() is None:
            log.warn("Failed to parse {}".format(mgafile['fn']))
        elif content.getroot().tag != 'MultiGenomeAlignmentSummaries':
            log.debug("{}/{} is not an MGA summary file.".format(mgafile['root'], mgafile['fn']))
        else:
            log.debug("Found file {}/{}".format(mgafile["root"], mgafile["fn"]))
            runId = content.findtext("RunID")
            if runId not in self.mga_data:
                self.mga_data[runId] = content


    def get_bar_data(self, mga_data):
        
        bar_data = OrderedDict()
        categories = OrderedDict()
        
        for summary in mga_data.findall("/MultiGenomeAlignmentSummary"):
            dataset_id = summary.findtext("DatasetId")
            sequence_count = int(summary.findtext("SequenceCount"))
            sampled_count = int(summary.findtext("SampledCount"))
            adapter_count = int(summary.findtext("AdapterCount"))
            unmapped_count = int(summary.findtext("UnmappedCount"))
            
            samples = self.get_sample_information(summary)
            
            species = set()
            controls = set()
            for sample in samples:
                if 'species' in sample:
                    sp = sample['species']
                    species.add(sp)
                    if sample['control']:
                        controls.add(sp)
            
            log.info("Dataset {} Species: {}".format(dataset_id, species))
            log.info("Dataset {} Controls: {}".format(dataset_id, controls))
            
            dataset_bar_data = dict()
            dataset_categories = dict()
            
            for alignmentSummary in summary.findall("AlignmentSummaries/AlignmentSummary"):
                assigned_count = int(alignmentSummary.findtext("AssignedCount"))
                
                assigned_error_rate = float(alignmentSummary.findtext("AssignedErrorRate"))
                error_rate = float(alignmentSummary.findtext("ErrorRate"))
                
                reference_genome_id = alignmentSummary.find("ReferenceGenome").attrib['id']
                reference_genome_name = alignmentSummary.find("ReferenceGenome").attrib['name']
                
                category_id = "{}.{}".format(dataset_id, reference_genome_id)
                
                colour = bar_colours.red
                if reference_genome_name in controls:
                    colour = bar_colours.orange
                elif reference_genome_name in species:
                    colour = bar_colours.green
                elif len(species) == 0 or 'other' in map(lambda s: s.lower(), species):
                    colour = bar_colours.grey
                
                alpha = max_alpha - (max_alpha - min_alpha) * (assigned_error_rate - min_error) / (max_error - min_error)
                alpha = max(min_alpha, min(max_alpha, alpha))
                
                # colour = colour.applyAlpha(alpha)
                
                if assigned_count >= 100:
                    log.debug("{}\t{}\t{}\t{}".format(reference_genome_id, assigned_count, error_rate * 100.0, alpha))
                    dataset_bar_data[category_id] = assigned_count
                
                    dataset_categories[category_id] = {
                        'name': reference_genome_name,
                        'color': colour.toHtml()
                    }
                
            # Sort into decreasing order of assigned count.
            dataset_bar_data = OrderedDict(sorted(dataset_bar_data.items(), key = lambda x: -x[1]))
            
            # The order of categories matters for the order in the plot, so add them to the
            # overall categories dictionary in the order they are in for the bar data.
            for category_id in dataset_bar_data.keys():
                categories[category_id] = dataset_categories[category_id]
                    
            log.info("Adapter count: {} / {}".format(adapter_count, sampled_count))
            
            dataset_bar_data['adapter'] = adapter_count

            bar_data[dataset_id] = dataset_bar_data

        categories['adapter'] = {
            'name': 'Adapter',
            'color': bar_colours.adapter.toHtml()
        }
        
        log.info("Bar data = {}".format(bar_data))
        log.info("Categories = {}".format(categories))

        return bar_data, categories

    
    def get_table_data(self, summary):
        headers = OrderedDict()
        '''
        headers['reference_id'] = {
            'title': 'Reference ID',
            'description': 'Internal reference genome identifier',
            'scale': False
        }
        '''
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

        table_config = {
            'namespace': 'mga',
            'id': 'mga_stats_table.{}'.format(dataset_id),
            'table_title': dataset_id,
            'col1_header': 'Reference ID',
            'no_beeswarm': True
        }
        
        table_data = dict()
        
        sequence_count = int(summary.findtext('SequenceCount'))
        sampled_count = int(summary.findtext('SampledCount'))
        adapter_count = int(summary.findtext('AdapterCount'))
        unmapped_count = int(summary.findtext('UnmappedCount'))
        
        for alignmentSummary in summary.findall("AlignmentSummaries/AlignmentSummary"):
            aligned_count = int(alignmentSummary.findtext("AlignedCount"))
            aligned_error = float(alignmentSummary.findtext("ErrorRate"))
            unique_count = int(alignmentSummary.findtext("UniquelyAlignedCount"))
            unique_error = float(alignmentSummary.findtext("UniquelyAlignedErrorRate"))
            preferred_count = int(alignmentSummary.findtext("PreferentiallyAlignedCount"))
            preferred_error = float(alignmentSummary.findtext("PreferentiallyAlignedErrorRate"))
            assigned_count = int(alignmentSummary.findtext("AssignedCount"))
            assigned_error = float(alignmentSummary.findtext("AssignedErrorRate"))
            
            reference_genome_id = alignmentSummary.find("ReferenceGenome").attrib['id']
            reference_genome_name = alignmentSummary.find("ReferenceGenome").attrib['name']
        
            aligned_fraction = float(aligned_count) / float(sampled_count)
            assigned_fraction = float(assigned_count) / float(sampled_count)
            
            # TODO. Got to here. More filters in results.xsl.
            
            if aligned_fraction >= aligned_fraction_threshold and assigned_fraction >= assigned_fraction_threshold:
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
        
        # TODO : Add others, controls, adapter
        
        return table.plot(table_data, headers, table_config)
        

    def get_sample_information(self, summary):
        samples = []
        for sample in summary.findall("Samples/Sample"):
            sample_info = { 'control': 'No' }
            self.read_properties(sample, sample_info)
            sample_info['control'] = sample_info['control'].lower() in ['yes','true']
            samples.append(sample_info)
        return samples

    
    def read_properties(self, element, props = dict()):
        for prop in element.findall("Properties/Property"):
            try:
                props[prop.attrib['name'].replace(' ', '').lower()] = prop.attrib['value']
            except KeyError:
                # Leave put any that have no value.
                pass
        return props


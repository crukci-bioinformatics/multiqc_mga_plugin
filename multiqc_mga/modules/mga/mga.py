#!/usr/bin/env python

""" MultiQC Multi Genome Alignment module """

import locale
import logging
import operator
import os

from collections import OrderedDict
from functools import cmp_to_key
from lxml import objectify, etree
from types import SimpleNamespace

from multiqc import config
from multiqc.plots import bargraph, table
from multiqc.modules.base_module import BaseMultiqcModule

from .colour import Colour

log = config.logger

# See https://stackoverflow.com/a/16279578
bar_colours = SimpleNamespace(**{
    'reference': Colour.fromBytes(0, 255, 0),
    'control': Colour.fromBytes(255, 200, 0),
    'contaminant': Colour.fromBytes(255, 0, 0),
    'unmapped': Colour.fromBytes(255, 255, 255),
    'unknown': Colour.fromBytes(128, 128, 128),
    'adapter': Colour.fromBytes(255, 102, 255)
})

max_alpha = 1.0
min_alpha = 0.4 # Is 0.1 in the Java version.
max_error = 0.01
min_error = 0.0025

assigned_fraction_threshold = 0.01
aligned_fraction_threshold = 0.01
error_rate_threshold = 0.0125

# Based on https://github.com/MultiQC/example-plugin

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name = 'Multi Genome Alignment',
            target = "MGA",
            anchor = 'mga',
            href = "https://github.com/crukci-bioinformatics/MGA",
            info = """(multi-genome alignment) is a quality control tool for high-throughput sequence data
                   written by Matthew Eldridge at the Cancer Research UK Cambridge Institute."""
        )

        # Add to self.css and self.js to be included in template
        self.css = { 'assets/css/multiqc_mga.css' : os.path.join(os.path.dirname(__file__), 'assets', 'css', 'multiqc_mga.css') }

        # XML files we know we don't want to load and check as they're from Illumina sequencing programs.
        # ConversionStats.xml in particular is quite large, so we don't want to parse it.
        self.illumina_files = [ 'ConversionStats.xml', 'DemultiplexingStats.xml', 'RunInfo.xml',
                                'runParameters.xml', 'RunParameters.xml', 'LaserPowerVariability.xml' ]

        self.sum_sequences = etree.XPath("sum(//MultiGenomeAlignmentSummaries/MultiGenomeAlignmentSummary/SequenceCount)")
        self.count_references = etree.XPath("count(//MultiGenomeAlignmentSummaries/ReferenceGenomes/ReferenceGenome)")

        mga_data = self.load_mga_files()
        self.create_mga_reports(mga_data)


    def load_mga_files(self):
        mga_data = dict()
        for mgafile in self.find_log_files('mga', filecontents=False, filehandles=True):
            with mgafile['f'] as fh:
                try:
                    if mgafile['fn'] not in self.illumina_files:
                        self._read_mga_xml_file(mga_data, mgafile)
                finally:
                    fh.close()

        if len(mga_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(mga_data)))

        return mga_data


    def _read_mga_xml_file(self, mga_data, mgafile):
        content = objectify.parse(mgafile["f"])
        if content is None or content.getroot() is None:
            log.warn("Failed to parse {}".format(mgafile['fn']))
        elif content.getroot().tag != 'MultiGenomeAlignmentSummaries':
            log.debug("{}/{} is not an MGA summary file.".format(mgafile['root'], mgafile['fn']))
        else:
            log.debug("Found file {}/{}".format(mgafile["root"], mgafile["fn"]))
            run_id = content.findtext("RunId")
            if run_id not in mga_data:
                mga_data[run_id] = content


    def create_mga_reports(self, mga_data):
        for run_id, dataset in mga_data.items():
            dataset_props = self._read_properties(dataset)

            total_sequence_count = self.sum_sequences(dataset)
            yield_multiplier = 2 if dataset_props.get('endtype') == 'Paired End' else 1
            cycles = int(dataset_props['cycles']) if 'cycles' in dataset_props else 0
            total_yield = yield_multiplier * cycles * total_sequence_count / 1000000000.0

            bar_data, bar_categories, plot_config = self.get_bar_data(dataset)

            bargraph_plot_html = bargraph.plot(bar_data, bar_categories, plot_config)

            trim_start = int(float(dataset.findtext("TrimStart")))
            trim_length = int(float(dataset.findtext("TrimLength")))
            number_of_genomes = int(self.count_references(dataset))

            self.add_section(
                name = 'Sequencing Results',
                description = self._plot_description_markdown(dataset, total_sequence_count, total_yield),
                content = self._plot_key_html(),
                anchor = 'mga_plot',
                helptext = '''
                    Sequences were sampled, trimmed to {} bases starting from position {}, and mapped to {} reference genomes
                    ([see list below](#mga_reference_genomes)) using Bowtie. Sequences containing adapters were found by
                    ungapped alignment of the full length sequence to a set of known adapter and primer sequences using
                    Exonerate. Further details on the alignment results and the assignment of reads to genomes are given
                    [below](#mga_alignment_details).
                '''.format(trim_length, trim_start, number_of_genomes),
                plot = bargraph_plot_html
            )

            sequencing_dataset = self._is_sequencing_dataset(dataset)

            for summary in dataset.findall("MultiGenomeAlignmentSummary"):
                dataset_id = summary.findtext("DatasetId")

                dataset_props = self._read_properties(dataset)

                sequence_count = int(summary.findtext('SequenceCount'))
                sampled_count = int(summary.findtext('SampledCount'))
                dataset_yield = yield_multiplier * cycles * float(sequence_count) / 1000000000.0

                self.add_section(
                    name = 'Lane {} Statistics'.format(dataset_id) if sequencing_dataset else 'Dataset "{}" Statistics'.format(dataset_id),
                    anchor = "mga_stats_{}".format(dataset_id),
                    plot = self.get_table_data(summary),
                    description = '''
                        <table>
                            <tr><td>Yield (Gbases):</td><td>{:.2f}</td></tr>
                            <tr><td>Sequences:</td><td>{:,}</td></tr>
                            <tr><td>Sampled:</td><td>{:,}</td></tr>
                        </table>
                        <br/>
                    '''.format(dataset_yield, sequence_count, sampled_count)
                )

                self.add_section(
                    name = 'Lane {} Samples'.format(dataset_id) if sequencing_dataset else 'Dataset "{}" Statistics'.format(dataset_id),
                    anchor = "mga_samples_{}".format(dataset_id),
                    plot = self.get_sample_table_data(summary)
                )

            # See https://stackoverflow.com/questions/2440692/formatting-floats-without-trailing-zeros
            self.add_section(
                name = "Alignment Details",
                anchor = "mga_alignment_details",
                description = """
                    Reference genomes are sorted according to how many sequence reads have been assigned to each.
                    Separate entries are given for reference genomes for which at least {:.4g}% of reads have been assigned
                    or for which at least {:.4g}% of reads align with an average mismatch or error rate of below {:.4g}%.

                    In addition to the total number of reads aligning to each reference genome and the average error
                    rate for those alignments, details are also provided for the the number of reads aligning uniquely
                    to the reference genome and and the associated error rate for those unique reads.

                    The 'Best' column and accompanying error rate refer to those reads that align preferentially to
                    the given reference genome, i.e. with the fewest mismatches. These reads will include those that
                    align uniquely and those that also align to other genomes with the same number of mismatches but
                    which do not align to another genome with fewer mismatches.

                    Reads that align uniquely to a genome are assigned to that genome. Reads that align equally well to
                    multiple genomes are assigned to the genome with the highest number of reads in the 'Best' column.

                    Note that because reads are trimmed prior to alignment with Bowtie, it is possible for a read to be
                    counted both as aligned to one or more of the reference genomes and among the reads with adapter
                    content. The adapter will most likely be present in the portion of the read that has been trimmed.
                """.format(assigned_fraction_threshold * 100, aligned_fraction_threshold * 100, error_rate_threshold * 100)
            )

            self.add_section(
                name = "Reference Genomes",
                anchor = "mga_reference_genomes",
                description = self.list_reference_genomes(dataset)
            )

    def get_bar_data(self, dataset):

        run_id = dataset.findtext("RunId")
        if run_id is None:
            raise UserWarning

        bar_data = OrderedDict()
        categories = OrderedDict()
        max_sequenced_count = 0

        for summary in dataset.findall("/MultiGenomeAlignmentSummary"):
            dataset_id = summary.findtext("DatasetId")
            sequence_count = int(summary.findtext("SequenceCount"))
            sampled_count = int(summary.findtext("SampledCount"))
            adapter_count = int(summary.findtext("AdapterCount"))
            unmapped_count = int(summary.findtext("UnmappedCount"))

            sampled_to_sequenced = float(sequence_count) / float(sampled_count)

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

                    colour = bar_colours.contaminant
                    if reference_genome_name in controls:
                        colour = bar_colours.control
                    elif reference_genome_name in species:
                        colour = bar_colours.reference
                    elif len(species) == 0 or 'other' in map(lambda s: s.lower(), species):
                        colour = bar_colours.unknown

                    # log.debug("Dataset {} - {}".format(dataset_id, reference_genome_id))
                    # log.debug("max_alpha - (max_alpha - min_alpha) * (assigned_error - min_error) / (max_error - min_error)")
                    # log.debug("{} - ({} - {}) * ({} - {}) / ({} - {})".format(max_alpha, max_alpha, min_alpha, assigned_error, min_error, max_error, min_error))
                    # log.debug("{} - {} * {} / {}".format(max_alpha, max_alpha - min_alpha, assigned_error - min_error, max_error - min_error))

                    alpha = max_alpha - (max_alpha - min_alpha) * (assigned_error - min_error) / (max_error - min_error)

                    # log.debug("alpha = {}".format(alpha))

                    alpha = max(min_alpha, min(max_alpha, alpha))

                    # log.debug("capped alpha = {}".format(alpha))

                    if assigned_count >= 100:
                        log.debug("{}\t{}\t{}\t{}".format(reference_genome_id, assigned_count, aligned_error * 100.0, alpha))

                    dataset_bar_data[category_id] = int(assigned_count * sampled_to_sequenced)

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

            dataset_bar_data['adapter'] = int(adapter_count * sampled_to_sequenced)

            bar_data[dataset_id] = dataset_bar_data

        categories['adapter'] = {
            'name': 'Adapter',
            'color': bar_colours.adapter.toHtml()
        }

        log.debug("Bar data = {}".format(bar_data))
        log.debug("Categories = {}".format(categories))

        plot_config = {
            'id': "mga_plot_{}".format(run_id.replace(' ', '_')),
            'title': "Multi Genome Alignment: {}".format(run_id),
            'cpswitch_counts_label': 'Number of reads',
            'xlab': "Lane" if self._is_sequencing_dataset(dataset) else "Data set",
            'ylab': "Number of reads",
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
            'min': 0,
            'max': 100,
            'format': '{:,.1%}',
            'scale': False
        }
        headers['aligned_error'] = {
            'title': 'Error rate',
            'description': 'Aligned error rate',
            'min': 0,
            'max': 100,
            'format': '{:,.2%}',
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
            'min': 0,
            'max': 100,
            'format': '{:,.2%}',
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
            'min': 0,
            'max': 100,
            'format': '{:,.2%}',
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
            'min': 0,
            'max': 100,
            'format': '{:,.2%}',
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
                    'aligned_perc': aligned_fraction,
                    'aligned_error': aligned_error,
                    'unique_count': unique_count,
                    'unique_error': unique_error,
                    'preferred_count': preferred_count,
                    'preferred_error': preferred_error,
                    'assigned_count': assigned_count,
                    'assigned_perc': assigned_fraction,
                    'assigned_error': assigned_error
                }

        # Sort into decreasing order of assigned count.
        table_data = OrderedDict(sorted(table_data.items(), key = lambda x: -x[1]['assigned_count']))

        if number_of_others >= 2:
            table_data['Other'] = {
                'species': "{} others".format(number_of_others),
                'aligned_count': other_assigned_count,
                'aligned_perc': other_assigned_count / sampled_count
            }

        table_data['Unmapped'] = {
            'species': '',
            'aligned_count': unmapped_count,
            'aligned_perc': unmapped_count / sampled_count
        }

        table_data['Adapter'] = {
            'species': '',
            'aligned_count': adapter_count,
            'aligned_perc': adapter_count / sampled_count
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

            # Make sure no element is None. Some keys may be in the sample properties
            # with no value, which gives None even if the default is used above.
            for key, value in table_data[row_id].items():
                if value is None:
                    table_data[row_id][key] = ''

        table_config = {
            'namespace': 'mga',
            'id': 'mga_sample_table.{}'.format(dataset_id),
            'table_title': dataset_id,
            'col1_header': 'Sample ID and name',
            'no_beeswarm': True,
            'sortRows': False
        }

        return table.plot(table_data, headers, table_config)


    def list_reference_genomes(self, mga_data):

        genomes = []
        for genome in mga_data.findall("ReferenceGenomes/ReferenceGenome"):
            genomes.append(genome.attrib['name'])

        # See https://stackoverflow.com/questions/36139/how-to-sort-a-list-of-strings
        genomes = sorted(genomes, key = cmp_to_key(locale.strcoll))

        markdown = "Sequences were aligned to the following reference genomes ({} in total).\n\n".format(len(genomes))

        for genome in genomes:
            markdown += "* {}\n".format(genome)

        return markdown


    def _is_sequencing_dataset(self, dataset):
        dataset_props = self._read_properties(dataset)
        return 'flowcellid' in dataset_props or 'runname' in dataset_props


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
                props[prop.attrib['name'].replace(' ', '').lower()] = prop.attrib.get('value')
            except KeyError:
                # Leave put any that have no name.
                pass
        return props


    def _plot_description_markdown(self, dataset, total_sequence_count, total_yield):
        results_desc = "<table>"
        for prop in dataset.findall("Properties/Property"):
            try:
                results_desc += "<tr><td>{}:</td><td>{}</td></tr>".format(prop.attrib['name'], prop.attrib.get('value'))
            except KeyError:
                pass
        results_desc += "<tr><td>Yield (Gbases):</td><td>{:.2f}</td></tr>".format(total_yield)
        results_desc += "<tr><td>Total sequences:</td><td>{:,.0f}</td></tr>".format(total_sequence_count)
        results_desc += "</table><br/>"

        return results_desc


    def _plot_key_html(self):
        key_elem_span = '<span class="multiqc_mga_key_element" style="background-color:{}">&nbsp;&nbsp;&nbsp;&nbsp;</span>&nbsp;{}\n'

        key_desc = '<div id="mga_plot_key">\n'
        key_desc += key_elem_span.format(bar_colours.reference.toHtml(), 'Sequenced&nbsp;species/genome')
        key_desc += key_elem_span.format(bar_colours.control.toHtml(), 'Control')
        key_desc += key_elem_span.format(bar_colours.contaminant.toHtml(), 'Contaminant')
        key_desc += key_elem_span.format(bar_colours.adapter.toHtml(), 'Adapter')
        key_desc += key_elem_span.format(bar_colours.unmapped.toHtml(), 'Unmapped')
        key_desc += key_elem_span.format(bar_colours.unknown.toHtml(), 'Unknown')
        key_desc += '</div>\n'

        return key_desc

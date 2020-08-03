#!/usr/bin/env python

""" MultiQC Multi Genome Alignment module """

import locale
import logging
import operator
import os

from collections import OrderedDict
from functools import cmp_to_key
from lxml import objectify, etree
from natsort import natsorted
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

gigabase = 1e9

assigned_fraction_threshold = 0.01
aligned_fraction_threshold = 0.01
error_rate_threshold = 0.0125
adapter_threshold_multiplier = 0.005

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

        # Add to self.css to be included in template
        self.css = { 'assets/css/multiqc_mga.css' : os.path.join(os.path.dirname(__file__), 'assets', 'css', 'multiqc_mga.css') }

        # XML files we know we don't want to load and check as they're from Illumina sequencing programs.
        # ConversionStats.xml in particular is quite large, so we don't want to parse it.
        self.illumina_files = [ 'ConversionStats.xml', 'DemultiplexingStats.xml', 'RunInfo.xml',
                                'runParameters.xml', 'RunParameters.xml', 'LaserPowerVariability.xml' ]

        mga_data = self.load_mga_files()

        self.create_mga_reports(mga_data)


    def load_mga_files(self):
        """Find and load the content of MGA XML files.

        :return A dictionary of run id (MGA, not necessary sequencing) to a list of XML
        trees with content for that run.
        :rtype dict
        """

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
        """Read a potential MGA XML file.

        This method has to load the content to make sure the file is from MGA. Ideally,
        MultiQC should be configured to have a more constrained name than "*.xml".
        Suitable content is added to the mga_data dictionary. Multiple content for the
        same run is added to the list stored in the dictionary.

        :param dict mga_data: A dictionary into which found content can be added, keyed by run id.
        :param dict mgafile: The XML file load which may or may not be an MGA file.
        """

        if mgafile['fn'] not in self.illumina_files:
            content = objectify.parse(mgafile["f"])
            if content is None or content.getroot() is None:
                log.warn("Failed to parse {}".format(mgafile['fn']))
            elif content.getroot().tag != 'MultiGenomeAlignmentSummaries':
                log.debug("{}/{} is not an MGA summary file.".format(mgafile['root'], mgafile['fn']))
            else:
                log.debug("Found file {}/{}".format(mgafile["root"], mgafile["fn"]))
                run_id = content.findtext("RunId")
                datasets = mga_data.get(run_id)
                if datasets is None:
                    datasets = []
                    mga_data[run_id] = datasets
                datasets.append(content)


    def _merge_run_summaries(self, run_datasets):
        '''
        Merge the XML trees for a single run id into a more useful form.
        
        :param run_datasets: The list of MGA summaries (XML trees) whose run ids are all the same.
        
        :return A dictionary that contains the summary of the information gathered for the run,
        all the MGA summary information into a dictionary keyed by dataset id,
        and all the reference genomes used.
        '''
        run_info = None
        mga_summaries_by_id = dict()
        reference_genomes = dict()

        for content in run_datasets:
            run_name = content.findtext("RunId")
            assert run_name is not None

            if run_info is None:
                summary_props = self._read_properties(content)

                run_info = {
                    'run_id': run_name,
                    'trim_start': int(float(content.findtext("TrimStart"))),
                    'trim_length': int(float(content.findtext("TrimLength"))),
                    'properties': summary_props,
                    'from_sequencing': 'Flow Cell ID' in summary_props or 'Run name' in summary_props,
                    'reference_genomes': reference_genomes,
                    'yield_multiplier': 2 if summary_props.get('End type') == 'Paired End' else 1,
                    'cycles': int(summary_props.get('Cycles', '0')),
                    'max_sequence_count': 0,
                    'total_sequence_count': 0
                }

            for mga_summary in content.findall("MultiGenomeAlignmentSummary"):
                dataset_id = mga_summary.findtext("DatasetId")
                assert dataset_id is not None
                if dataset_id in mga_summaries_by_id:
                    log.warn(f'MGA run "{run_name}" has more than one occurrence of mga_summary "{dataset_id}".')
                else:
                    mga_summaries_by_id[dataset_id] = mga_summary

                    sequence_count = int(mga_summary.findtext("SequenceCount"))
                    run_info['max_sequence_count'] = max(run_info['max_sequence_count'], sequence_count)
                    run_info['total_sequence_count'] = run_info['total_sequence_count'] + sequence_count

            for reference in content.findall("ReferenceGenomes/ReferenceGenome"):
                reference_genomes[reference.attrib['id']] = reference.attrib['name']

        run_info['total_yield'] = run_info['yield_multiplier'] * run_info['cycles'] * run_info['total_sequence_count'] / gigabase

        # Sort the datasets into a natural order for their dataset id.
        run_info['mga_summaries'] = OrderedDict(natsorted(mga_summaries_by_id.items(), key = lambda x: x[0]))

        run_info = SimpleNamespace(**run_info)

        log.debug("Run {} has {} summaries and {} reference genomes.".format(
                  run_info.run_id, len(run_info.mga_summaries), len(run_info.reference_genomes)))

        return run_info


    def create_mga_reports(self, mga_data):
        for run_id, run_contents in mga_data.items():
            run_info = self._merge_run_summaries(run_contents)

            plot_data, plot_categories = self._plot_data(run_info)

            self.add_section(
                name = 'Sequencing Results',
                anchor = 'mga_sequencing_results',
                description = self._plot_description(run_info),
                content = self._plot_key(),
                helptext = self._plot_help(run_info),
                plot = bargraph.plot(plot_data, plot_categories, self._plot_config(run_info))
            )

            for dataset_id, mga_summary in run_info.mga_summaries.items():

                sequence_count = int(mga_summary.findtext('SequenceCount'))
                sampled_count = int(mga_summary.findtext('SampledCount'))
                dataset_yield = run_info.yield_multiplier * run_info.cycles * float(sequence_count) / gigabase

                self.add_section(
                    name = f'Lane {dataset_id} Statistics' if run_info.from_sequencing else f'Dataset "{mga_summary}" Statistics',
                    anchor = f"mga_stats_{dataset_id}",
                    plot = table.plot(self._main_table_data(mga_summary), self._main_table_headers(), self._main_table_config(dataset_id)),
                    description = f"""
                        <table>
                            <tr><td>Yield (Gbases):</td><td>{dataset_yield:.2f}</td></tr>
                            <tr><td>Sequences:</td><td>{sequence_count:,}</td></tr>
                            <tr><td>Sampled:</td><td>{sampled_count:,}</td></tr>
                        </table>
                        <br/>
                    """
                )

                self.add_section(
                    name = f'Lane {dataset_id} Samples' if run_info.from_sequencing else f'Dataset "{dataset_id}" Statistics',
                    anchor = f"mga_samples_{dataset_id}",
                    plot = table.plot(self._sample_table_data(mga_summary), self._sample_table_headers(), self._sample_table_config(run_info, dataset_id))
                )


    def _plot_data(self, run_info):

        bar_data = OrderedDict()
        categories = OrderedDict()

        for dataset_id, mga_summary in run_info.mga_summaries.items():
            sequence_count = int(mga_summary.findtext("SequenceCount"))
            sampled_count = int(mga_summary.findtext("SampledCount"))
            adapter_count = int(mga_summary.findtext("AdapterCount"))
            unmapped_count = int(mga_summary.findtext("UnmappedCount"))

            sampled_to_sequenced = float(sequence_count) / float(sampled_count)

            species, controls = self._get_species_and_controls(mga_summary)

            log.debug("Dataset {} Species: {}".format(dataset_id, species))
            log.debug("Dataset {} Controls: {}".format(dataset_id, controls))

            dataset_bar_data = dict()
            dataset_categories = dict()

            for alignment_summary in mga_summary.findall("AlignmentSummaries/AlignmentSummary"):
                if self._accept_genome(species, mga_summary, alignment_summary):
                    aligned_count = int(alignment_summary.findtext("AlignedCount"))
                    aligned_error = float(alignment_summary.findtext("ErrorRate"))
                    assigned_count = int(alignment_summary.findtext("AssignedCount"))
                    assigned_error = float(alignment_summary.findtext("AssignedErrorRate"))

                    reference_genome_id = alignment_summary.find("ReferenceGenome").attrib['id']
                    reference_genome_name = alignment_summary.find("ReferenceGenome").attrib['name']

                    category_id = f"{dataset_id}.{reference_genome_id}"

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

            bar_data[dataset_id] = dataset_bar_data

            # The order of categories matters for the order in the plot, so add them to the
            # overall categories dictionary in the order they are in for the bar data.
            for category_id in dataset_bar_data.keys():
                categories[category_id] = dataset_categories[category_id]

            log.debug(f"Adapter count: {adapter_count} / {sampled_count}")

            if adapter_count >= sampled_count * adapter_threshold_multiplier:
                dataset_adapter_id = f"{dataset_id}A" if run_info.from_sequencing else f"{dataset_id} adapter"
                category_id = f"{dataset_id}.adapter"
                dataset_bar_data = { category_id: int(adapter_count * sampled_to_sequenced) }
                bar_data[dataset_adapter_id] = dataset_bar_data
                categories[category_id] = {
                    'name': 'Adapter',
                    'color': bar_colours.adapter.toHtml()
                }


        log.debug(f"Bar data = {bar_data}")
        log.debug(f"Categories = {categories}")

        return bar_data, categories


    def _plot_config(self, run_info):
        return {
            'id': "mga_plot_{}".format(run_info.run_id.replace(' ', '_')),
            'title': f"Multi Genome Alignment: {run_info.run_id}",
            'cpswitch_counts_label': 'Number of reads',
            'xlab': "Lane" if run_info.from_sequencing else "Data set",
            'ylab': "Number of reads",
            'ymin': 0,
            'ymax': run_info.max_sequence_count,
            'use_legend': False,
            'tt_percentages': False
        }


    def _plot_description(self, run_info):
        results_desc = "<table>"
        for name, value in run_info.properties.items():
            try:
                results_desc += "<tr><td>{}:</td><td>{}</td></tr>".format(name, value)
            except KeyError:
                pass
        results_desc += f"<tr><td>Yield (Gbases):</td><td>{run_info.total_yield:.2f}</td></tr>"
        results_desc += f"<tr><td>Total sequences:</td><td>{run_info.total_sequence_count:,.0f}</td></tr>"
        results_desc += "</table><br/>"

        return results_desc


    def _plot_key(self):
        key_elem_span = '<span class="multiqc_mga_key_element" style="background-color:{}">&nbsp;&nbsp;&nbsp;&nbsp;</span>&nbsp;{}\n'

        key_desc = '<div id="mga_plot_key">\n'
        key_desc += key_elem_span.format(bar_colours.reference.toHtml(), 'Sequenced&nbsp;species/genome')
        key_desc += key_elem_span.format(bar_colours.control.toHtml(), 'Control')
        key_desc += key_elem_span.format(bar_colours.contaminant.toHtml(), 'Contaminant')
        key_desc += key_elem_span.format(bar_colours.unmapped.toHtml(), 'Unmapped')
        key_desc += key_elem_span.format(bar_colours.unknown.toHtml(), 'Unknown')
        key_desc += key_elem_span.format(bar_colours.adapter.toHtml(), 'Adapter')
        key_desc += '</div>\n'

        return key_desc


    def _plot_help(self, run_info):

        # See https://stackoverflow.com/questions/36139/how-to-sort-a-list-of-strings
        genomes = natsorted(run_info.reference_genomes.values(), key = cmp_to_key(locale.strcoll))
        number_of_genomes = len(genomes)

        help = f"""
            Sequences were sampled, trimmed to {run_info.trim_length} bases starting from position {run_info.trim_start},
            and mapped to {number_of_genomes} reference genomes (see list below) using Bowtie. Sequences containing
            adapters were found by ungapped alignment of the full length sequence to a set of known adapter and
            primer sequences using Exonerate.

        """

        if run_info.from_sequencing:
            help += f"""
                Some lanes may have an addition bar with the suffix "A" (e.g. "1A" for lane 1).
                This is the adapter when the number of adapter reads is significant, equal to or
                above {adapter_threshold_multiplier * 100:.4g}% of the number of sampled reads.
            """
        else:
            help += f"""
                Some data sets may have an addition bar with the suffix "adapter".
                This is the adapter when the number of adapter reads is significant for the
                data set, equal to or above {adapter_threshold_multiplier * 100:.4g}% of the number of sampled reads.
            """

        # See https://stackoverflow.com/questions/2440692/formatting-floats-without-trailing-zeros

        help += f"""

            #### Alignment Details

            Reference genomes are sorted according to how many sequence reads have been assigned to each.
            Separate entries are given for reference genomes for which at least {assigned_fraction_threshold * 100:.4g}%
            of reads have been assigned or for which at least {aligned_fraction_threshold * 100:.4g}% of reads align
            with an average mismatch or error rate of below {error_rate_threshold * 100:.4g}%.

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

            #### Reference Genomes

            Sequences were aligned to the following reference genomes ({number_of_genomes} in total).

        """

        for genome in genomes:
            help += f"* {genome}\n"

        return self._strip_from_lines(help)


    def _main_table_data(self, mga_summary):

        dataset_id = mga_summary.findtext('DatasetId')
        sequence_count = int(mga_summary.findtext('SequenceCount'))
        sampled_count = int(mga_summary.findtext('SampledCount'))
        adapter_count = int(mga_summary.findtext('AdapterCount'))
        unmapped_count = int(mga_summary.findtext('UnmappedCount'))

        species, controls = self._get_species_and_controls(mga_summary)

        number_of_others = 0
        other_assigned_count = 0
        table_data = dict()

        for alignment_summary in mga_summary.findall("AlignmentSummaries/AlignmentSummary"):
            if not self._accept_genome(species, mga_summary, alignment_summary):
                assigned_count = int(alignment_summary.findtext("AssignedCount"))
                number_of_others = number_of_others + 1
                other_assigned_count = other_assigned_count + assigned_count

        for alignment_summary in mga_summary.findall("AlignmentSummaries/AlignmentSummary"):
            if number_of_others < 2 or self._accept_genome(species, mga_summary, alignment_summary):
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
                'species': f"{number_of_others} others",
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

        return table_data


    def _main_table_headers(self):
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
            'scale': False,
            #'shared_key': 'read_count'
        }
        headers['aligned_perc'] = {
            'title': 'Aligned %',
            'description': 'Percentage of reads aligned',
            'min': 0,
            'max': 100,
            'format': '{:,.1%}',
            'scale': False,
            #'shared_key': 'percent_aligned'
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
            'scale': False,
            #'shared_key': 'read_count'
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
            'scale': False,
            #'shared_key': 'read_count'
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
            'scale': False,
            #'shared_key': 'read_count'
        }
        headers['assigned_perc'] = {
            'title': 'Assigned %',
            'description': 'Percentage of reads assigned',
            'suffix': '%',
            'min': 0,
            'max': 100,
            'format': '{:,.1f}',
            'scale': False,
            #'shared_key': 'percent_aligned'
        }
        headers['assigned_error'] = {
            'title': 'Error rate',
            'description': 'Assigned error rate',
            'min': 0,
            'max': 100,
            'format': '{:,.2%}',
            'scale': False
        }
        return headers


    def _main_table_config(self, dataset_id):
        return {
            'namespace': 'mga',
            'id': f'mga_stats_table.{dataset_id}',
            'table_title': dataset_id,
            'col1_header': 'Reference ID',
            'no_beeswarm': True,
            'sortRows': False
        }


    def _sample_table_data(self, mga_summary):
        dataset_id = mga_summary.findtext('DatasetId')

        table_data = OrderedDict()

        for sample_xml in mga_summary.findall("Samples/Sample"):
            sample_id = None
            sample_name = None
            sample_info = dict()
            for name, value in self._read_properties(sample_xml).items():
                key = name.lower().replace(' ', '_')
                if key == 'sample_id':
                    sample_id = value
                elif key == 'sample_name':
                    sample_name = value
                else:
                    # Make sure no element is None. Some keys may be in the sample properties
                    # with no value.
                    sample_info[key] = '' if value is None else value

            row_id = f"{sample_id} / {sample_name}"
            table_data[row_id] = sample_info

        return table_data


    def _sample_table_headers(self):
        headers = OrderedDict()
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
        return headers


    def _sample_table_config(self, run_info, dataset_id):
        return {
            'namespace': 'mga',
            'id': f'mga_sample_table.{dataset_id}',
            'table_title': dataset_id,
            'col1_header': "Pool / Sample" if run_info.from_sequencing else "Sample ID / Name",
            'no_beeswarm': True,
            'sortRows': False,
            'only_defined_headers': True
        }


    def _get_species_and_controls(self, summary):
        samples = self._get_sample_information(summary)

        species = set()
        controls = set()
        for sample in samples:
            if 'Species' in sample:
                sp = sample['Species']
                species.add(sp)
                if sample.get('Control', False):
                    controls.add(sp)

        return species, controls


    def _accept_genome(self, species_set, mga_summary, alignment_summary):
        reference_genome_name = alignment_summary.find("ReferenceGenome").attrib['name']

        if reference_genome_name in species_set:
            return True

        sampled_count = int(mga_summary.findtext("SampledCount"))

        aligned_count = int(alignment_summary.findtext("AlignedCount"))
        aligned_error_rate = float(alignment_summary.findtext("ErrorRate"))
        assigned_count = int(alignment_summary.findtext("AssignedCount"))

        aligned_fraction = float(aligned_count) / float(sampled_count)
        assigned_fraction = float(assigned_count) / float(sampled_count)

        return assigned_fraction >= assigned_fraction_threshold or aligned_fraction >= aligned_fraction_threshold and aligned_error_rate < error_rate_threshold


    def _get_sample_information(self, mga_summary):
        samples = []
        for sample in mga_summary.findall("Samples/Sample"):
            sample_info = { 'Control': 'No' }
            self._read_properties(sample, sample_info)
            sample_info['Control'] = sample_info['Control'].lower() in ['yes','true']
            samples.append(sample_info)
        return samples


    def _read_properties(self, element, props = OrderedDict()):
        for prop in element.findall("Properties/Property"):
            try:
                props[prop.attrib['name']] = prop.attrib.get('value')
            except KeyError:
                # Leave out any that have no name.
                pass
        return props

    def _strip_from_lines(self, str):
        lines = [ l.strip() for l in str.splitlines() ]
        return "\n".join(lines)


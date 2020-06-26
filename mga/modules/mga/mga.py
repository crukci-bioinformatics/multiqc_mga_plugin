#!/usr/bin/env python

""" MultiQC Multi Genome Alignment module """

import logging
import operator
import os
import untangle
from collections import OrderedDict, defaultdict
from itertools import islice
from lxml import objectify
from types import SimpleNamespace

from multiqc import config
from multiqc.plots import bargraph, table
from multiqc.modules.base_module import BaseMultiqcModule

log = config.logger

# See https://stackoverflow.com/a/16279578
bar_colours = SimpleNamespace(**{
    'red': '#ff0000',
    'orange': '#ffcb00',
    'green': '#00ff00',
    'grey': '#808080',
    'adapter': '#ff66ff'
})

max_alpha = 1.0
min_alpha = 0.1
max_error = 0.0025
min_error = 0.01

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
        
        for run in self.mga_data.values():
            for summary in run.findall("/MultiGenomeAlignmentSummary"):
                self.get_bar_data_from_summary(summary)

    def read_mga_xml_file(self, mgafile):
        #content = untangle.parse(mgafile["f"])
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

    def get_bar_data_from_summary(self, summary):
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
        
        log.debug("Species: {}".format(species))
        log.debug("Controls: {}".format(controls))
        
        bar_data = dict()
        
        for alignmentSummary in summary.findall("AlignmentSummaries/AlignmentSummary"):
            assigned_count = int(alignmentSummary.findtext("AssignedCount"))
            
            assigned_error_rate = float(alignmentSummary.findtext("AssignedErrorRate"))
            error_rate = float(alignmentSummary.findtext("ErrorRate"))
            
            reference_genome_id = alignmentSummary.find("ReferenceGenome").attrib['id']
            reference_genome_name = alignmentSummary.find("ReferenceGenome").attrib['name']
            
            colour = bar_colours.red
            if reference_genome_name in controls:
                colour = bar_colours.orange
            elif reference_genome_name in species:
                colour = bar_colours.green
            elif len(species) == 0 or 'other' in map(lambda s: s.lower(), species):
                colour = bar_colours.grey
            
            alpha = max_alpha - (max_alpha - min_alpha) * (assigned_error_rate - min_error) / (max_error - min_error)
            alpha = max(min_alpha, min(max_alpha, alpha))
            
            if assigned_count >= 100:
                log.debug("{}\t{}\t{}\t{}".format(reference_genome_id, assigned_count, error_rate * 100.0, alpha))
            
            if assigned_count > 0:
                bar_data[reference_genome_id] = {
                    'genome': reference_genome_name,
                    'count': assigned_count,
                    'colour': colour
                }
            
        # Sort into decreasing order of assigned count.
        bar_data = OrderedDict(sorted(bar_data.items(), key = lambda x: -x[1]['count']))
                
        log.debug("Adapter count: {} / {}".format(adapter_count, sampled_count))
        
        bar_data['adapter'] = {
            'count': adapter_count,
            'colour': bar_colours.adapter
        }
        
        log.info("Bar data = {}".format(bar_data))

        return bar_data

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
            props[prop.attrib['name'].replace(' ', '').lower()] = prop.attrib['value']
        return props


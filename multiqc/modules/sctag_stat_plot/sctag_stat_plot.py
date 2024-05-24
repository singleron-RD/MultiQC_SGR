""" MultiQC module to parse output from tsne """

import logging
import os
import csv
from collections import defaultdict

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Stat",
            anchor="stat"
        )

        # Find and load any sctagPlot reports
        self.sctag_data = defaultdict(lambda:defaultdict(dict))
        self.plot_data = defaultdict(lambda:defaultdict(dict))
        ## map
        for f in self.find_log_files("sctag_stat_plot/map",filehandles=True):
            self.parse_mapplot_report(f,f["f"],f["s_name"])
        ## count
        for f in self.find_log_files("sctag_stat_plot/count",filehandles=True):
            self.parse_countplot_report(f,f["f"],f["s_name"])

        # Filter to strip out ignored sample names
        self.sctag_data = self.ignore_samples(self.sctag_data)

        if len(self.sctag_data) == 0:
            raise ModuleNoSamplesFound
        num_rep = 0
        for _step in self.sctag_data:
            for _sample in self.sctag_data[_step]:
                num_rep += 1
        log.info(f"Found {num_rep} reports")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Write parsed report data to a file
        for _step in self.sctag_data:
            self.write_data_file(self.sctag_data[_step], f"multiqc_{_step}")

        # Assignment bar plot
        self.add_section(
            name = 'Mapping',
            anchor = 'map',
            helptext = """Reads Mapped : R2 reads that successfully mapped to linker and tag-barcode.
Reads Unmapped Invalid Linker : Unmapped R2 reads because of too many mismatches in linker sequence.
Reads Unmapped Invalid Barcode : Unmapped R2 reads because of too many mismatches in tag-barcode sequence.""",
            plot=self.plot_chart('map'))
        self.add_section(
            name = 'Cells',
            anchor = 'count',
            plot=self.plot_chart('count'))
    
    def parse_mapplot_report(self,f,ff,s_name):
        """Parse the Plot log file."""

        parsed_data = dict()
        s_name = s_name.removesuffix("_stat")

        reader = csv.reader(ff,delimiter=":")
        for row in reader:
            parsed_data[row[0].strip()] = row[1].strip()

        data = dict()     
        for k in parsed_data:
            data[k] = int(parsed_data[k].split("(")[0].replace(',',''))
                
        self.plot_data['map'] = data
        # Add to the main dictionary
        if len(data) > 1:
            if s_name in self.sctag_data['map']:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            self.add_data_source(f, s_name)
            self.sctag_data['map'][s_name] = data

    def parse_countplot_report(self,f,ff,s_name):
        """Parse the Plot log file."""

        parsed_data = dict()
        s_name = s_name.removesuffix("_stat")

        reader = csv.reader(ff,delimiter=":")
        for row in reader:
            parsed_data[row[0].strip()] = row[1].strip()

        data = dict()     
        for k in parsed_data:
            if k not in ('Number of Matched Cells','Mapped Reads in Cells','Median UMI per Cell','Mean UMI per Cell'):
                data[k] = int(parsed_data[k].split("(")[0].replace(',',''))
                
        self.plot_data['count'] = data
        # Add to the main dictionary
        if len(data) > 1:
            if s_name in self.sctag_data['count']:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            self.add_data_source(f, s_name)
            self.sctag_data['count'][s_name] = data

    
    def plot_chart(self,step):
        """Make the Plot assignment rates plot"""

        headers = {}
        for h in self.plot_data[step]:
            headers[h] = {"name": h}

        # Config for the plot
        config = None
        if step == 'map':
            config = {
                "id": "extract_plot",
                "title": "Mapping sequence proportion",
                "ylab": "# Reads",
                "cpswitch_counts_label": "Number of Reads",
            }
        elif step == 'count':
            config = {
                "id": "extract_plot",
                "title": "The proportion of tags in cells",
                "ylab": "# Cells",
                "cpswitch_counts_label": "Number of Cells",
            }

        return bargraph.plot(self.sctag_data[step], headers, config)
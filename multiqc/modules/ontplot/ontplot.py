""" MultiQC module to parse output from ontPlot """


import logging
import re
import csv

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="ONT extract",
            anchor="ontplot",
            info='Statistical chart of ONT sequence splitting.'
        )

        # Find and load any ontPlot reports
        self.ontplot_data = dict()
        self.plot_data = dict()
        for f in self.find_log_files("ontplot/plot",filehandles=True):
            self.parse_ontplot_report(f,f["f"],f["s_name"])

        # Filter to strip out ignored sample names
        self.ontplot_data = self.ignore_samples(self.ontplot_data)

        if len(self.ontplot_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.ontplot_data)} reports")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Write parsed report data to a file
        self.write_data_file(self.ontplot_data, "multiqc_ontPlot")

        # Assignment bar plot
        self.add_section(plot=self.ontplot_chart())

    def parse_ontplot_report(self,f,ff,s_name):
        """Parse the ontPlot log file."""

        parsed_data = dict()
        s_name = s_name.removesuffix(".Summary")

        reader = csv.reader(ff)
        for row in reader:
            parsed_data[row[0]] = row[1]
    
        data = dict()     
        for k in parsed_data:
        
            if (k.find('valid') != -1) | (k.find('frac') != -1):
                continue
            else:
                data[k] = parsed_data[k]
        self.plot_data = data
        del self.plot_data['total']
        # Add to the main dictionary
        if len(data) > 1:
            if s_name in self.ontplot_data:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            self.add_data_source(f, s_name)
            self.ontplot_data[s_name] = data

    
    def ontplot_chart(self):
        """Make the ontPlot assignment rates plot"""

        headers = {}
        for h in self.plot_data:
            headers[h] = {"name": h}

        # Config for the plot
        config = {
            "id": "ont_extract_plot",
            "title": "Split sequence proportion",
            "ylab": "# Reads",
            "cpswitch_counts_label": "Number of Reads",
        }
        return bargraph.plot(self.ontplot_data, headers, config)

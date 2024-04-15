""" MultiQC module to parse output from STARSolo """


import logging
import json

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph


# Initialise the logger
log = logging.getLogger(__name__)


def get_frac(x) -> str:
    x = round(float(x) * 100, 2)
    return str(x)


def get_int(x) -> str:
    return str(int(x))


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="STARSolo",
            anchor="starsolo",
            href="https://github.com/alexdobin/STAR",
            info="for scRNA-Seq",
            doi="10.1093/bioinformatics/bts635",
        )

        # Find and load any STAR reports
        summary_data = self.parse_log("summary")
        read_stats_data = self.parse_log("read_stats")
        if len(summary_data) == 0 and len(read_stats_data) == 0:
            raise ModuleNoSamplesFound

        # Basic Stats Table
        self.general_stats_table(summary_data)
        # plot
        self.add_section(
            name="Read Counts Assigned to Features", anchor="starsolo_assign", plot=self.assign_chart(read_stats_data)
        )

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

    def parse_log(self, seg):
        data_dict = dict()
        for f in self.find_log_files(f"starsolo/{seg}"):
            parsed_data = json.loads(f["f"])
            if parsed_data is not None:
                s_name = f["s_name"].removesuffix(f".{seg}")
                if s_name in data_dict:
                    log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                self.add_data_source(f, s_name=s_name, section=seg)
                data_dict[s_name] = parsed_data

        data_dict = self.ignore_samples(data_dict)

        log.info(f"Found {len(data_dict)} starsolo {seg} reports")
        # Write parsed report data to a file
        self.write_data_file(data_dict, f"multiqc_starsolo_{seg}")
        return data_dict

    def general_stats_table(self, summary_data):
        headers = {
            "Reads With Valid Barcodes": {
                "title": "% Valid Barcodes",
                "description": "Fraction of reads with valid barcodes matching whitelist",
                "format": get_frac,
                "max": 100,
                "min": 0,
                "suffix": "%",
            },
            "Estimated Number of Cells": {
                "title": "N Cells",
                "description": "Estimated number of cells",
                "format": get_int,
            },
            "Fraction of Unique Reads in Cells": {
                "title": "% Reads in Cells",
                "description": "Fraction of unique reads in cells",
                "format": get_frac,
                "max": 100,
                "min": 0,
                "suffix": "%",
            },
            "Median GeneFull_Ex50pAS per Cell": {
                "title": "Median Genes",
                "description": "Median genes per cell",
                "format": get_int,
            },
            "Mean Reads per Cell": {
                "title": "Mean Reads",
                "description": "Mean Reads per Cell",
                "format": get_int,
            },
            "Mean UMI per Cell": {
                "title": "Mean UMI",
                "description": "Mean UMI per Cell",
                "format": get_int,
            },
            "Sequencing Saturation": {
                "title": "Saturation",
                "description": "Sequencing Saturation",
                "format": get_frac,
                "max": 100,
                "min": 0,
                "suffix": "%",
            },
        }
        self.general_stats_addcols(summary_data, headers=headers)

    def assign_chart(self, read_stats_data):
        """Make the plot showing alignment rates"""

        # Specify the order of the different possible categories
        keys = {
            "exonic": {"color": "#437bb1", "name": "Mapped reads assigned to exonic regions"},
            "intronic": {"color": "#7cb5ec", "name": "Mapped reads assigned to intronic regions"},
            "intergenic": {"color": "#e63491", "name": "Mapped reads assigned to intergenic regions"},
            "antisense": {"color": "#f7a35c", "name": "Mapped reads assigned antisense to gene"},
        }

        # Config for the plot
        pconfig = {
            "id": "starsolo_assign_plot",
            "title": "STARSolo: Assign",
            "ylab": "# Reads",
            "cpswitch_counts_label": "Number of Reads",
        }

        return bargraph.plot(read_stats_data, keys, pconfig)

""" MultiQC module to parse output from ONT_extract """


import logging
import os
import csv

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound


# Initialise the logger
log = logging.getLogger(__name__)

def get_frac(x) -> str:
    x = round(float(x) * 100, 2)
    return str(x) + "%"

def get_int(x) -> str:
    return str(int(x))


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="ONT_stat",
            anchor="ontstat"
        )

        # Find and load any STAR reports
        self.ontextract_data = dict()
        for f in self.find_log_files("ontstat", filehandles=True):
            parsed_data = self.parse_ontextract_report(f["f"])
            if parsed_data is not None:
                s_name = f["s_name"]
                s_name = s_name.removesuffix(".csv")
                if s_name in self.ontextract_data:
                    log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                self.add_data_source(s_name=s_name, source=os.path.join(f["root"], f["fn"]), section="summary")
                self.ontextract_data[s_name] = parsed_data

        self.ontextract_data = self.ignore_samples(self.ontextract_data)
        if len(self.ontextract_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.ontextract_data)} reports")
        # Write parsed report data to a file
        self.write_data_file(self.ontextract_data, "multiqc_ontstat")

        # Basic Stats Table
        self.general_stats_table()

    def parse_ontextract_report(self, f):
        parsed_data = {}
        reader = csv.reader(f)
        for row in reader:
            parsed_data[row[0]] = row[1]
        return parsed_data

    def general_stats_table(self):
        headers = {
            "Estimated Number of Cells": {
                "title": "Estimated Number of Cells",
                "description": "The number of barcodes associated with at least one cell.",
                "format": get_int,
            },
            "Mean Reads per Cell": {
                "title": "Mean Reads per Cell",
                "description": "The total number of sequenced reads divided by the number of barcodes associated with cell-containing partitions.",
                "format": get_int,
            },
            "Median Genes per Cell": {
                "title": "Median Genes per Cell",
                "description": "The median number of genes detected per cell-associated barcode. Detection is defined as the presence of at least 1 UMI count.",
                "format": get_int,
            },
            "% Valid Barcodes": {
                "title": "% Valid Barcodes",
                "description": "Fraction of reads with valid barcodes matching whlitelist.",
                "format": get_frac,
            },
            "% Reads in Cells": {
                "title": "% Reads in Cells",
                "description": "Fraction of unique reads in cells",
                "format": get_frac,
            },
            "Mean UMI per Cell": {
                "title": "Mean UMI per Cell",
                "description": "Mean UMI per Cell",
                "format": get_int,
            },

        }
        self.general_stats_addcols(self.ontextract_data, headers=headers)

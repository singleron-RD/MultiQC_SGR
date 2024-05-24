""" MultiQC module to parse output from tsne """

import logging
import os
import csv
from collections import defaultdict

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound


# Initialise the logger
log = logging.getLogger(__name__)

def get_frac(x) -> str:
    x = round(float(x), 2)
    return str(x)

def get_int(x) -> str:
    return str(int(x))

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Stat",
            anchor="stat"
        )

        # Find and load any STAR reports
        self.extract_data = defaultdict(dict)
        num_rep = 0
        for f in self.find_log_files("sctag_stat_table/bc",filehandles=True):
            num_rep += 1
            parsed_data = self.parse_bcextract_report(f["f"])
            if parsed_data is not None:
                s_name = f["s_name"].removesuffix("_stat")
                self.add_data_source(s_name=s_name, source=os.path.join(f["root"], f["fn"]), section="summary")
                self.extract_data[s_name].update(parsed_data)
        for f in self.find_log_files("sctag_stat_table/cutadapt",filehandles=True):
            num_rep += 1
            parsed_data = self.parse_cutextract_report(f["f"])
            if parsed_data is not None:
                s_name = f["s_name"].removesuffix("_stat")
                self.add_data_source(s_name=s_name, source=os.path.join(f["root"], f["fn"]), section="summary")
                self.extract_data[s_name].update(parsed_data)

        self.extract_data = self.ignore_samples(self.extract_data)
        if len(self.extract_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {num_rep} reports")
        # Write parsed report data to a file
        self.write_data_file(self.extract_data, "multiqc_bc_cutadapt")

        # Basic Stats Table
        self.general_stats_table()

    def parse_bcextract_report(self, f):
        parsed_data = {}
        reader = csv.reader(f,delimiter=":")
        for row in reader:
            _use_name = row[0].strip()
            _use_data = row[1].strip()
            if _use_name in ('Raw Reads','Valid Reads','Q30 of Barcodes','Q30 of UMIs'):
                if _use_data.find("(") != -1:
                    _use_data = _use_data.split('(')[1].split('%')[0]
                else:
                    _use_data = _use_data.split('%')[0].replace(',','')
                parsed_data[_use_name] = _use_data
        return parsed_data
    
    def parse_cutextract_report(self, f):
        parsed_data = {}
        reader = csv.reader(f,delimiter=":")
        for row in reader:
            _use_name = row[0].strip()
            _use_data = row[1].strip().split('(')[1].split('%')[0]
            parsed_data[_use_name] = _use_data
        return parsed_data

    def general_stats_table(self):
        headers = {
            "Raw Reads": {
                "title": "Raw Reads",
                "description": "Total reads from FASTQ files.",
                "format": get_int,
            },
            "Valid Reads": {
                "title": "% Valid Reads",
                "description": "Reads pass filtering(filtered: reads without poly T, reads without linker, reads without correct barcode or low quality reads).",
                "format": get_frac,
                'max': 100,
                'min': 0,
                'scale': 'OrRd',
                'suffix': '%'
            },
            "Q30 of Barcodes": {
                "title": "% Q30 of Barcodes",
                "description": "Percent of barcode base pairs with quality scores over Q30.",
                "format": get_frac,
                'max': 100,
                'min': 0,
                'scale': 'RdPu',
                'suffix': '%'
            },
            "Q30 of UMIs": {
                "title": "% Q30 of UMIs",
                "description": "Percent of UMI base pairs with quality scores over Q30.",
                "format": get_frac,
                'max': 100,
                'min': 0,
                'scale': 'YlOrRd',
                'suffix': '%'
            },
            "Reads with Adapters": {
                "title": "% Reads with Adapters",
                "description": "Reads with sequencing adapters or reads two with poly A(read-through adpaters).",
                "format": get_frac,
                'max': 100,
                'min': 0,
                'scale': 'YlOrRd-rev',
                'suffix': '%'
            },
            "Reads Too Short": {
                "title": "% Reads Too Short",
                "description": "Reads with read length less than 20bp after trimming.",
                "format": get_frac,
                'max': 100,
                'min': 0,
                'scale': 'RdPu-rev',
                'suffix': '%'
            },
            "Base Pairs Quality-Trimmed": {
                "title": "% Base Pairs Quality-Trimmed",
                "description": "Bases pairs removed from the end of the read whose quality is smaller than the given threshold.",
                "format": get_frac,
                'max': 100,
                'min': 0,
                'scale': 'OrRd-rev',
                'suffix': '%'
            }
        }
        self.general_stats_addcols(self.extract_data, headers=headers)
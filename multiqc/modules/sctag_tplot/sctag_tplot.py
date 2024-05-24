""" MultiQC module to parse output from tsne """

import logging
import json

from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="t-SNE plot",
            anchor="tsne_plot"
        )
        self.tsne_cluster_data = dict()
        self.tsne_gene_data = dict()
        self.s_name_set = set()
        # Find and load any tsne reports
        for i,f in enumerate(self.find_log_files("sctag_tplot",filehandles=True),1):
            self.parse_tsne_report(f["f"],f["s_name"])

        log.info(f"Found {i} reports")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Assignment plot
        for sample in self.s_name_set:
            self.add_section(
                plot=self.tsne_cluster_data[sample]
            )
        for sample in self.s_name_set:
            self.add_section(
                plot=self.tsne_gene_data[sample]
            )

    def parse_tsne_report(self,ff,s_name):
        """Parse the tsne log file."""
        _name = s_name.removesuffix("_html")
        _data = json.load(ff)
        self.s_name_set.add(_name)
        self.tsne_cluster_data[_name] = _data['tsne_cluster']
        self.tsne_gene_data[_name] = _data['tsne_gene']

""" MultiQC module to parse output from plbcrank """


import logging
import csv

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import linegraph
from multiqc.modules.cellranger.utils import parse_bcknee_data

import plotly.graph_objs as go
import plotly.offline as pltoff
import math


# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Barcode rank",
            anchor="plbcrank"
        )

        # Find and load any plbcrank reports
        self.plbcrank_conf = {"bc": dict()}
        self.plbcrank_data = {"bc": dict()}
        self.plot_data = dict()
        self.chart = None
        for f in self.find_log_files("plbcrank",filehandles=True):
            self.parse_plbcrank_report(f,f["f"],f["s_name"])

        # Filter to strip out ignored sample names
        for k in self.plbcrank_data.keys():
            self.plbcrank_data[k] = self.ignore_samples(self.plbcrank_data[k])

        log.info(f"Found {len(self.plbcrank_data)} reports")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Write parsed report data to a file
        self.write_data_file(self.plbcrank_data, "multiqc_plbcrank")

        # Assignment bar plot
        self.add_section(
            name="Count - BC rank plot",
            anchor="count-bcrank-plot",
            description=self.plbcrank_conf["bc"]["description"],
            helptext=self.plbcrank_conf["bc"]["helptext"],
            plot=self.chart
        )


    def parse_plbcrank_report(self,f,ff,s_name):
        """Parse the plbcrank log file."""
        x_list = []
        y_list = []
        s_name = s_name.removesuffix(".bcrank")
        
        reader = csv.reader(ff)
        next(reader)
        for i,row in enumerate(reader,1):
            x_list.append(i)
            y_list.append(int(row[1]))

        anchor_cell_umi = int(y_list[30]/10)
        y_anchor = len([y for y in y_list if y > anchor_cell_umi])

        cell_num = y_anchor
        plot_data_list = [{'x': x_list[:y_anchor],
                            'y': y_list[:y_anchor],
                            'type': 'scattergl',
                            'mode': 'lines',
                            'line': {'color': '#4682B4', 'width': 3},
                            'showlegend': True,
                            'name': 'Cells',
                            'text': f'100% Cells<br>({cell_num}/{cell_num})'},
                          {'x': x_list[y_anchor-1:],
                            'y': y_list[y_anchor-1:],
                            'type': 'scattergl',
                            'mode': 'lines',
                            'line': {'color': '#dddddd', 'width': 3},
                            'showlegend': True,
                            'name': 'Background',
                            'text': f'Background<br>'}
                            ]
        layout = go.Layout(width=950, height=650,
                title={"text": "Barcode rank", "font": {"color": "black"}, "x": 0.5},
                xaxis={"type": "log", "title": "Barcode", "titlefont": {"color": "black"},
                        "color": "black", "gridcolor": "gainsboro", "linecolor": "black"},
                yaxis={"type": "log", "title": "UMI counts", "titlefont": {"color": "black"},
                        "color": "black", "gridcolor": "gainsboro", "linecolor": "black"},
                margin=dict(l=50, r=0, t=30, b=30),
                plot_bgcolor="#FFFFFF")
        config = dict({"displayModeBar": True,
                    "staticPlot": False,
                    "showAxisDragHandles": False,
                    "modeBarButtons": [["toImage", "resetScale2d"]],
                    "scrollZoom": False,
                    "displaylogo": False})

        fig = go.Figure(data=plot_data_list, layout=layout)
        self.chart = pltoff.plot(fig, include_plotlyjs=True, output_type='div', config=config)

        plots = {
            "bc": {
                "config": {
                    "id": "mqc_cellranger_count_bc_knee",
                    "title": "Barcode Rank Plot",
                    "xlab": 'Barcodes',
                    "ylab": 'UMI counts',
                    "yLog": True,
                    "xLog": True,
                },
                "description": "Barcode knee plot",
                "helptext": 'The plot shows the count of filtered UMIs mapped to each barcode.  As barcodes are not determined to be cell-associated strictly based on their UMI count, but instead are determined by their expression profiles, some regions of the graph contain both cell-associated and background-associated barcodes.  The color of the graph in these regions is based on the local density of barcodes that are cell-associated.',
            }}
        self.plbcrank_conf.update(plots)
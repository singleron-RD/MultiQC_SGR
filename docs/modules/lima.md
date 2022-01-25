---
Name: Lima
URL: https://github.com/PacificBiosciences/barcoding
Description: >
  Lima, the PacBio barcode demultiplexer, is the standard tool to
  identify barcode sequences in PacBio single-molecule sequencing data. Starting
  in SMRT Link v5.1.0, it is the tool that powers the Demultiplex Barcodes
  GUI-based analysis application.
---

The Lima module parses the report and count files generated by
[Lima](https://github.com/PacificBiosciences/barcoding), a PacBio tool to
demultiplex PacBio single-molecule sequencing data.

By default, Lima will use `barcode1--barcode2` as the sample names. To prevent
these barcodes showing up in the General Statistics table of MultiQC, the Lima
results are added to their own section.

If you want to include the Lima results in the General Statistics table, you
can rename the `barcode1--barcode2` filenames to their apropriate samples using
the [--replace-names](https://multiqc.info/docs/#sample-name-replacement)
option. Each sample that is specified in this way will be moved from the Lima
section to the General Statistics table.
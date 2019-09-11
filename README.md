# Telomere_sequence_content_analysis_tools
Collection of scripts developed to analyse telomere sequence content

## Scripts for telomere controls (telo_ctrl/)
### sort_telo_ctrls.pl
Usage: perl sort_telo_ctrls.pl <input.fastq>

Description: Sorts reads into 4 files, (1) reads that are from the G-rich telomere strand, (2) reads that are from the C-rich telomere strand, (3) reads that contain a mix of both G- and C-rich telomere repeats, and (4) any reads remaining reads that do not fit into the former 3 groups.

Output files: *.G-strands.fastq, *.C-strands.fastq, *.mixed-strands.fastq, *.others.fastq

### variant_analysis.pl
Usage: perl sort_telo_ctrls.pl <input.fastq>

Description: Analyses the reads from the fastq file, counting the number of each single nucleotide variant in the TTAGGG repeat unit and outputing the counts in csv format. The counts are given in long format, with the counts categorised by base type, base position in TTAGGG repeat unit, position in read, and the total length of the read.

Output files: *.variants_analysis_6.csv


# Metagenomic analysis of animal dental calculus
Details on code and bioinformatic tools used to analyse shotgun DNA sequencing datasets.
All the analyses were run in the HPC server Galielo (Cineca, Italy).

## Download of raw sequencing data
The raw sequencing data of the dataset analysed were downloaded with SRA Toolkit by generating a text file `SRA_Acc_List.txt` with the list of SRA accession number and the following commands: 

```bash
cat SRA_Acc_List.txt | xargs -I{} prefetch {}
fastq-dump --split-files --gzip *.sra
vdb-validate filenames.sra
```
Data from Weyrich et al. were downloaded from the Online Ancient Gene Repository (https://www.oagr.org.au/) with `wget`.

## Quality-filtering and adapter-trimming of the reads
All the dental calculus data consisted of paired-end sequencing reads, which were quality-filtered, trimmed of the adapter sequences and merged with AdapterRemoval: 

`AdapterRemoval --file1 filename_1.fastq.gz --file2 filename_2.fastq.gz --basename filename --minlength 30 --minquality 25 --trimns --trimqualities --gzip --threads 16 --collapse`




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

## Pre-processing of raw sequencing data
All the dental calculus data consisted of paired-end sequencing reads, which were quality-filtered, trimmed of the adapter sequences and merged with AdapterRemoval: 

```bash
AdapterRemoval --file1 filename_1.fastq.gz --file2 filename_2.fastq.gz --basename filename --minlength 30 --minquality 25 --trimns --trimqualities --gzip --threads 16 --collapse
```

We removed exact duplicates and reverse exact duplicates (`-derep 14`) from the merged reads with Prinseq: 
```bash
gzip -dc filename.fasta.gz | prinseq-lite.pl -verbose -fastq stdin -graph_data filename.prinseq.gd -graph_stats ld,de -derep 14 -out_good stdout -out_bad null | gzip > filename_prinseq.fastq.gz
```

## Taxonomic classification of the reads

### Construction of Kraken2 custom database
We constructed a custom database of 50 Gb including complete genomes of Bacteria, Archaea and Viruses from the NCBI RefSeq (as of November 2020) following the manual of Kraken2: 

```bash
DBNAME=customkraken2_50Gb_Nov2020
DBSIZE=50000000000
THREADS=16
kraken2-build --download-taxonomy --threads $THREADS --db $DBNAME
kraken2-build --download-library bacteria --threads $THREADS --db $DBNAME
kraken2-build --download-library viral --threads $THREADS --db $DBNAME
kraken2-build --download-library archaea --threads $THREADS --db $DBNAME
kraken2-build --build --threads $THREADS --db $DBNAME --max-db-size $DBSIZE 
```

### Taxonomic classification of reads with Kraken2 custom database
The custom database was used to make the taxonomic classification of the pre-processed reads. 
```bash
kraken2 --db $DBNAME --threads $THREADS filename.collapsed.gz --output filename.krk --gzip-compressed --report filename.krk.report
```

We ran Bracken to obtain actual quantification data, as described in the manual. We first built the Bracken database, and then ran the samples. 
We used a read length of 65 bp. 
```bash
KRAKEN_DB=customkraken2_50Gb_Nov2020
THREADS=16
READ_LEN=65
bracken -d ${KRAKEN_DB} -i filename.krk.report -o filename.bracken -r ${READ_LEN} 
```

We used the custom R script `brackenToAbundanceTable.R` to parse the Bracken abundance data for each individual in a table. 
```bash
DIR=path/to/bracken/output
brackenToAbundanceTable.R $DIR
```
The script generates two abundance tables, one with taxa as species names (as in the NCBI) and one with the species reported as NCBI IDs. For downastream analysis we used the table with species names. 





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

We used the custom R script `brackenToAbundanceTable.R` (see Toolbox repository) to parse the Bracken abundance data for each individual from a study dataset in a table. 
```bash
DIR=path/to/bracken/output
brackenToAbundanceTable.R $DIR
```
The script generates two abundance tables, one with taxa as species names (as in the NCBI) and one with the species reported as NCBI IDs. For downastream analysis we used the table with species names. 

### Filtration of contaminating sequences
In the baboon dental calculus dataset, the ancient Egyptian baboons were filtered of the contaminant sequences by selecting with `awk` in the abundance table the species with >10 and >200 reads, respectively, in the columns corresponding to the NTCs and the environmental controls (see details in Ottoni et al. 2019, Scientific Reports). 

```bash
awk -F'\t' -v OFS='\t' '($4 > 200) || ($7 > 200) || ($10 > 10) || ($11 > 10)  || ($12 > 10) || ($13 > 10){print $1}' taxa_abundance_bracken_names.txt | tail -n +2 > blanks10.env200.contaminants
```

A list of soil and skin species was generated from the abundance tables parsed from Bracken (present in two different folders, same filename), and they were merged in one file: 
```bash
cut -f1 taxa_abundance_bracken_names.txt > soil_species.list
cut -f1 taxa_abundance_bracken_names.txt > skin_species.list
cat skin_species.list soil_species.list | sort | uniq | grep -v "taxon" > soil_skin_species.contaminants
```

We downloaded the list of species in the HOMD, and used this list in a text file (`HOMD_species_Sept2020.txt`)to removed the species of possible oral origin present from the contaminants list:
```bash
grep -vFf HOMD_species_Sept2020.txt soil_skin_species.contaminants > soil_skin_species.homd.contaminants   
```

We merged the contaminant files and removed them from the abundance table: 
```bash
cat blanks10.env200.contaminants soil_skin_species.homd.contaminants | sort | uniq > full.contaminants
grep -vFf full.contaminants taxa_abundance_bracken_names > taxa_abundance_bracken_names.filtered
```

The same procedure was used to filter the data from Brealey et al. (2020, MBE). 

### Normalization for genome length and final abundance table
The abundance table of each dataset was normalized for genome lentgh with the custom Python script (see Toolbox repository). All the tables were merged in a single one with the custom R script `abundanceTablesMerger.R`.

```bash
gL-normalizer-lite.py taxa_abundance_bracken_names.txt prokaryotes_viruses_organelles.table taxa_abundance_bracken_names_norm.txt
abundanceTablesMerger.R *norm.txt
```

We generated an abundance table that included the full taxonomic classification (up to Phylum) for each species, which was then used to generate abudance table at the genus-level. We did that by using the script `getFullTaxaranks.sh` (see Toolbox repository). To run the script, you have to install the program taxonomy-ranks (https://github.com/linzhi2013/taxonomy_ranks/blob/master/README.md).

```bash
getFullTaxaranks.sh -i abundance_table.merged 
```

We removed the Virues from the table with `grep`. 

```bash
grep -v "Viruses" abundance_table.merged.final > abundance_table.merged.final.noVirus
```

After that, we generated a table of genus abundances with the script `genustabGenerator.R` (see Toolbox repository). The script creates a table `abundance_table.genus`

```bash
genustabGenerator.R abundance_table.merged.taxonomy.final.noVirus 
```

## R session
We imported the table in R and normalized the genus abundances for sequencing depth by reating a function for the Total Sum Scaling in which the read count of each genus is divided by the total number of read counts of the sample.

```R
tab = t(read.delim("abundance_table.genus", header=T, fill=T, row.names=1, sep="\t"))
TSS.divide = function(x){
 x/sum(x)
}
tab.tss = t(apply(tab.red, 1, TSS.divide))
```
We created a function to filter out the taxa represented below the threshold of 0.02% 

```R
low.count.removal = function(
                        data, # OTU count data frame of size n (sample) x p (OTU)
                        percent=0.02 # cutoff chosen
                        ){
    keep.otu = which(colSums(data)*100/(sum(colSums(data))) > percent)
    data.filter = data[,keep.otu]
    return(list(data.filter = data.filter, keep.otu = keep.otu))
}

result.filter = low.count.removal(tab.tss, percent=0.02)
# get the actual filtered data from results
tab.tss.flt = result.filter$data.filter
```

We created a metadata vector `group` that assigned to each sample in the columns of the table of genera abundances (normalized and filtered) the biological orgin of the sample, namely: Ancient baboon, Baboon zoo, Gorilla, Reindeer, Brown Bear, Chimpanzee, Historic Human, and so on (as defined in figure legends of the manuscript). 

### non-Metric Multidimensional Scaling (nMDS)

We used metaMDS in vegan to run the nMDS

```R
library(vegan)
nmds = metaMDS(tab.tss.flt, distance="bray", k=2, trymax=200, autotransform = FALSE, engine = "monoMDS", maxit = 200)
nmds = metaMDS(tab.tss.flt, distance="bray", k=2, trymax=200, autotransform = FALSE, engine = "monoMDS", maxit = 200, previous.best=nmds)
```

We set up `pch`, `bg` and `coul` to customize the charts and used plot and the metadata in the `group` vector to plot the results of the the nMDS

```R
plot(nmds, type="n", main="Bray-Curtis", cex.axis=0.75, cex.lab=0.75, xlab="", ylab="", xlim=c(-2.3,2.5), ylim=c(-2,2), xaxs="i", yaxs="i")
abline(v=0, lty=3, col="grey")
abline(h=0, lty=3, col="grey")
points(nmds, cex=1, bg=bg[factor(group)], pch=pch[factor(group)], col=coul[factor(group)], lwd=0.9)
```

### Cluster analysis (UPGMA)

We used vegan to calculated the Bray-Curtis and run the UPGMA with hclust

```R
library(vegan)
library(ape)
bray_dist = vegdist(tab.tss.flt, method = "bray", binary = FALSE, diag = FALSE, upper = FALSE, na.rm = FALSE)
bray_dist.clust = hclust(bray_dist, method="average", members = NULL)
plot(as.phylo(bray_dist.clust), type = "unrooted", cex = 0.6, lab4ut="axial", no.margin=T, show.tip.label=T, label.offset=0.02, edge.color = "gray", edge.width = 1, edge.lty = 1)
tiplabels(pch=pch[as.factor(group)], col = coul[as.factor(group)], cex=1.3, lwd=1, bg=bg[as.factor(group)])          
legend("topleft", legend = sort(unique(group)), bty = "n", col = coul, pch = pch, pt.cex=1.3, cex=0.8, pt.bg=bg, pt.lwd=1)
add.scale.bar(cex=0.8)
```

### Deseq2 analysis
We restricted the analysis to selected samples (numbering follows the table that we generated), converted the table in counts per million (cpm) and transformed abundances to integers (otherwise Deseq2 returns an error message). 

```R
tab.deseq2 = tab.tss.flt[c(1:6,9,10,		#Baboons
							13,16,			              #Gorillas
							17:21,			              #Reindeer
							22,23,24,26,	            #Bear
							28:46,255,		            #Chimps
							182:225,		              #Historic humans
							230:239),]		            #Modern humans
tab.deseq2.cpm = tab.deseq2*1000000
tab.deseq2.cpm.int <- tab.deseq2.cpm
for (i in 1:ncol(tab.deseq2.cpm)) {
    tab.deseq2.cpm.int[,i] <- as.integer(tab.deseq2.cpm[,i])
}
```

We created a metadata tab-separated text file in the form: 
Sample_ID\tGroup 

We imported in R the metadata file and run Deseq2 as follows: 

```R
metadata = read.delim("metadata_deseq2_animals.txt", header=T)

# move rownames to column 1 for deSeq format table.
library(DESeq2)
library(tibble)
tab.deseq2.cpm.int.ds = rownames_to_column(as.data.frame(t(tab.deseq2.cpm.int)), var = "species")

# Prepare deSeq
dds.data <- DESeqDataSetFromMatrix(countData=tab.deseq2.cpm.int.ds, 
                              colData=metadata, 
                              design=~group, tidy = TRUE)

# Run deSeq                              
dds.data = DESeq(dds.data)
```

Finally, we ran the contrast function of Deseq2 for each groups pair, we filtered for padj values and sorted the discriminant taxa based on the baseMean values.

``` R
res1 <- results(dds.data, contrast=c("group","Baboon","Reindeer"))
res2 <- results(dds.data, contrast=c("group","Baboon","Bear"))
res3 <- results(dds.data, contrast=c("group","Baboon","Chimpanzee"))
res4 <- results(dds.data, contrast=c("group","Baboon","Modern human"))	
res5 <- results(dds.data, contrast=c("group","Bear","Reindeer"))	
res6 <- results(dds.data, contrast=c("group","Bear","Chimpanzee"))	
res7 <- results(dds.data, contrast=c("group","Bear","Modern human"))
res8 <- results(dds.data, contrast=c("group","Reindeer","Chimpanzee"))
res9 <- results(dds.data, contrast=c("group","Reindeer","Modern human"))
res10 <- results(dds.data, contrast=c("group","Chimpanzee","Modern human"))
res11 <- results(dds.data, contrast=c("group","Historic human","Baboon"))
res12 <- results(dds.data, contrast=c("group","Historic human","Reindeer"))
res13 <- results(dds.data, contrast=c("group","Historic human","Bear"))
res14 <- results(dds.data, contrast=c("group","Historic human","Chimpanzee"))

# sort based on p-value adjusted
resOrdered1 <- res1[order(res1$padj),]
resOrdered2 <- res2[order(res2$padj),]
resOrdered3 <- res3[order(res3$padj),]
resOrdered4 <- res4[order(res4$padj),]
resOrdered5 <- res5[order(res5$padj),]
resOrdered6 <- res6[order(res6$padj),]
resOrdered7 <- res7[order(res7$padj),]
resOrdered8 <- res8[order(res8$padj),]
resOrdered9 <- res9[order(res9$padj),]
resOrdered10 <- res10[order(res10$padj),]
resOrdered11 <- res11[order(res11$padj),]
resOrdered12 <- res12[order(res12$padj),]
resOrdered13 <- res13[order(res13$padj),]
resOrdered14 <- res14[order(res14$padj),]

# get only significant taxa based on p-value adjusted
resSig1 <- subset(resOrdered1, padj < 0.1)
resSig2 <- subset(resOrdered2, padj < 0.1)
resSig3 <- subset(resOrdered3, padj < 0.1)
resSig4 <- subset(resOrdered4, padj < 0.1)
resSig5 <- subset(resOrdered5, padj < 0.1)
resSig6 <- subset(resOrdered6, padj < 0.1)
resSig7 <- subset(resOrdered7, padj < 0.1)
resSig8 <- subset(resOrdered8, padj < 0.1)
resSig9 <- subset(resOrdered9, padj < 0.1)
resSig10 <- subset(resOrdered10, padj < 0.1)
resSig11 <- subset(resOrdered11, padj < 0.1)
resSig12 <- subset(resOrdered12, padj < 0.1)
resSig13 <- subset(resOrdered13, padj < 0.1)
resSig14 <- subset(resOrdered14, padj < 0.1)

# sort significant values based on abundance:
x1 = resSig1[order(resSig1$baseMean),]
x2 = resSig2[order(resSig2$baseMean),]
x3 = resSig3[order(resSig3$baseMean),]
x4 = resSig4[order(resSig4$baseMean),]
x5 = resSig5[order(resSig5$baseMean),]
x6 = resSig6[order(resSig6$baseMean),]
x7 = resSig7[order(resSig7$baseMean),]
x8 = resSig8[order(resSig8$baseMean),]
x9 = resSig9[order(resSig9$baseMean),]
x10 = resSig10[order(resSig10$baseMean),]
x11 = resSig11[order(resSig11$baseMean),]
x12 = resSig12[order(resSig12$baseMean),]
x13 = resSig13[order(resSig13$baseMean),]
x14 = resSig14[order(resSig14$baseMean),]
```


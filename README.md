**Neranjan Perera**
<h1>get_orthogroupgenes</h1>

The program will selects the grouped ortholog genes for a particular species.		

usage: `get_orthogroup_genes.py [-h] [--orthofile ORTHO_GROUP_FILE] [-fasta FASTA_FILE] [--tag SPECIES_TAG] [-l] [-s]
                                [--outfile OUTPUT_FILE] [--version]`

optional arguments:
 - -h, --help            show this help message and exit
 
optional arguments:
 - -h, --help            show this help message and exit
 
optional arguments:
- -h, --help            show this help message and exit


usage: `get_orthogroup_genes.py [-h] [--orthofile ORTHO_GROUP_FILE]		
                               [--fasta FASTA_FILE] [--tag SPECIES_TAG] [-l]. 
                               [-s] [--outfile OUTPUT_FILE] [--version] `
                               

optional arguments: 
 ... -h, --help            show this help message and exit.
  
  --orthofile ORTHO_GROUP_FILE. 
                        ortho group file name. 
                        
  --fasta FASTA_FILE    FASTA file name. 
  
  --tag SPECIES_TAG     gene identification tag for the species eg: TR is the
                        tag for TR22118|c0_g1_i1_14120|m.21703. 
                        
  -l                    boolean switch selects the longest gene length. If not
                        used default is set to select the longest. 
                        
  -s                    boolean switch selects the shortest gene length. If
                        not used default is set to select the longest. 
                        
  --outfile OUTPUT_FILE
                        output file name of the sorted genes. If not given it
                        will use the default name; sorted_genes.fasta. 
                        
  --version             show program's version number and exit. 

s = "Python syntax highlighting"
print s.

--

```python
import sys,re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import MutableSeq
from array import array



def get_sequence(f_file,gid,out_file):
        f = open(f_file, 'r')
        for seq_record in SeqIO.parse(f,"fasta"):
                if ( gid == seq_record.id):
                        rec = SeqRecord(seq_record.seq,
                                id = seq_record.id,
                                description = "")
                        SeqIO.write(rec, out_file, "fasta")
                        break
        f.close()
        return;


def get_csv_geneID(csv_line):
        ln = csv_line.strip().split(",")
        return ln[1].split("\"")[1].split(":")[1]

def get_seqid(g_id,fasta_out):
        gtf = open(sys.argv[2],"r")
        for gtf_line in gtf:
                gtf_line = gtf_line.split('#', 1)[0]
                if(gtf_line != ""):
                        gln = gtf_line.strip().split("\t")
                        gid = gln[5]
                        if(gid == g_id):
                                print("GeneID:%s\ttableID:\t%s\tprotein-ID:\t%s" %(g_id, gid, gln[7]))
                                get_sequence(sys.argv[3],gln[7],fasta_out)
                                break

```

```sh
#!/bin/bash
#SBATCH --job-name=entap
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=50G
#SBATCH --partition=cbcworkshop
#SBATCH --qos=cbcworkshop
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo `hostname`

module load anaconda/2.4.0
module load perl/5.24.0
module load diamond/0.9.19
module load eggnog-mapper/0.99.1
module load interproscan/5.25-64.0

/UCHC/LABS/Wegrzyn/EnTAP/EnTAP --runP -i fasta_out.fasta -d /isg/shared/databases/Diamond/RefSeq/vertebrate_other.protein.faa.88.dmnd -d /isg/shared/databases/Diamond/Uniprot/uniprot_sprot.dmnd --ontology 0  --threads 16
```

```R
# Load DESeq2 library
library("DESeq2")

# Set the working directory
directory <- "~/Documents/R/DESeq2/"
setwd(directory)

# Set the prefix for each output file name
outputPrefix <- "Listeria_DESeq2"

# Location of deseq ready count table (EDGE-pro output)
deseqFile <- "~/Documents/R/DESeq2/Listeria_deseqFile"

# Read the table
countData <- read.table(deseqFile, header = T)

# Replace accession numbers with meaningful names
names(countData) <- c("10403S Rep1","DsigB Rep1","10403S Rep2","DsigB Rep2")

# Create table with treatment information
sampleNames <- colnames(countData)
sampleCondition <- c("10403S","DsigB","10403S","DsigB")
colData <- data.frame(condition = sampleCondition)
row.names(colData) = sampleNames
treatments = c("10403S","DsigB")

# Create DESeqDataSet: countData is the count table, colData is the table with treatment information
# One experimental condition
ddsFromMatrix <- DESeqDataSetFromMatrix(countData = countData,
 colData = colData,
 design = ~ condition)
colData(ddsFromMatrix)$condition <- factor(colData(ddsFromMatrix)$condition, levels = treatments)
dds <- DESeq(ddsFromMatrix)
res <- results(dds)

# filter results by p value
# order results by padj value (most significant to least)
res= subset(res, padj<0.05)
res <- res[order(res$padj),]
# should see DataFrame of baseMean, log2Foldchange, stat, pval, padj

# save data results and normalized reads to csv
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata)[1] <- 'gene'
write.csv(resdata, file = paste0(outputPrefix, "-results-with-normalized.csv"))

# send normalized counts to tab delimited file for GSEA, etc.
write.table(as.data.frame(counts(dds),normalized=T), file = paste0(outputPrefix, "_normalized_counts.txt"), sep = '\t')

# produce DataFrame of results of statistical tests
mcols(res, use.names = T)
write.csv(as.data.frame(mcols(res, use.name = T)),file = paste0(outputPrefix, "-test-conditions.csv"))

# replacing outlier value with estimated value as predicted by distrubution using
# "trimmed mean" approach. recommended if you have several replicates per treatment
# DESeq2 will automatically do this if you have 7 or more replicates
ddsClean <- replaceOutliersWithTrimmedMean(dds)
ddsClean <- DESeq(ddsClean)
tab <- table(initial = results(dds)$padj < 0.05,
 cleaned = results(ddsClean)$padj < 0.05)
addmargins(tab)
write.csv(as.data.frame(tab),file = paste0(outputPrefix, "-replaceoutliers.csv"))
resClean <- results(ddsClean)

# filter results by p value
resClean = subset(res, padj<0.05)
resClean <- resClean[order(resClean$padj),]
write.csv(as.data.frame(resClean),file = paste0(outputPrefix, "-replaceoutliers-results.csv"))
####################################################################################
# Exploratory data analysis of RNAseq data with DESeq2
#
# these next R scripts are for a variety of visualization, QC and other plots to
# get a sense of what the RNAseq data looks like based on DESEq2 analysis
#
# 1) MA plot
# 2) rlog stabilization and variance stabiliazation
# 3) heatmap of clustering analysis
# 4) PCA plot
#
#
####################################################################################

# MA plot of RNAseq data for entire dataset
# http://en.wikipedia.org/wiki/MA_plot
# genes with padj < 0.1 are colored Red
plotMA(dds, ylim=c(-8,8),main = "RNAseq Listeria monocytogenes")
dev.copy(png, paste0(outputPrefix, "-MAplot_initial_analysis.png"))
dev.off()
```

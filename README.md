# Approaches for Comparative Genomics of Replicate Microbial Communities

## Introduction

The aim of this workshop is to present a walk-through of the steps for
comparative genomic analysis of microbial communities. Most microbiome
studies use 16S rRNA sequencing to determine community composition. However,
any analysis based on conserved marker genes will by design only be able to
reach a limited level of genetic resolution such as genera or species. Most
microbes---bacteria in particular---also display a large amount of genetic
diversity with genera or species. This diversity is only accessible through
sequencing of whole genomes.

We have recently completed [a study](https://dx.doi.org/10.1101/280271) in
which we used whole-genome untargeted "shotgun" metagenomic sequencing to
describe the genetic diversity across a large number of replicate microbial
communities. The study focused on millimeter-scale granular microbial biofilms
that make up granular activated sludge.

## Full Analysis Overview

From samples to results follow these steps:

1. Sample preparation
    1. Extract DNA from samples
    2. Prepare libraries and send to sequencing

2. Quality filtering
    1. Pairing of paired-end reads
    2. Trimming of unpaired reads

3. Reference database
    - Unsupervised: _de-novo_ assembly of metagenome assembled genomes
      (MAGs)
    - Supervised: published genomes

4. Align short sequencing reads to reference database
    1. Construct database for aligner from multi FASTA file
    2. Run aligner (Bowtie2, BWA, MagicBLAST).  
       **Warning: this step is resource intensive and might be run on a compute cluster**

5. Derive community composition (i.e. "OTU table")
    1. Group reference genomes based on criteria (species, genus, family, etc.)
    2. Correct counts for genome length
    3. **Construct relative abundance table**
    4. Perform common analysis (e.g. in R)

6. [Genomic analysis of a target organism group](#markdown-header-6-genomic-analysis-of-accumulibacter)
    1. Extract reads that map to any of the genomes in the group
    2. Get reprensentative genomes
    3. Map reads to each of the genomes
    4. Get diversity at each base position
    5. Analyze variability within and between genomes
    6. Build consensus genomes
    7. Reconstruct phylogenies

## Required software

- Some Unix, e.g. Linux, macOS, or Cygwin
- [Bowtie2](http://bowtie-bio.sf.net/bowtie2)  
- [samtools](http://www.htslib.org)
- [bcftools](http://www.htslib.org)

These packages are really easy to install with common package managers.

- macOS
    - Homebrew: `brew install bowtie2 samtools bcftools`
- Linux
    - Linuxbrew: `brew install bowtie2 samtools bcftools`
    - Ubuntu: `apt-get install bowtie2 samtools bcftools`
    - CentOS: `yum install bowtie2 samtools bcftools`
- and more...

Additionally, we will require a couple of tools that are not released
in a formal package. These tools can be found as part of the
[gentools](https://www.bitbucket.org/gaberoo/gentools) repository.

```bash

git clone https://www.bitbucket.org/gaberoo/gentools.git
cd gentools

# for macOS
cp Make.inc.macos Make.inc

# for linux
cp Make.inc.linux Make.inc

make vcfConsensus vcf2hs vcf2fst

```

## 6  Genomic Analysis of Accumulibacter

For this workshop we will focus on the workflow for the building community-level consensus genomes, and then performing comparative analysis across them.


### 6.1  Extract reads that map to Accumulibacter

_Note, that the FASTQ files can be [downloaded
directly](https://www.dropbox.com/sh/k9biu9qladix882/AAARsIZNxf4bl6yk3LUUBnSxa
?dl=0) as well._

We assume that we have already successfully mapped all the short reads to the
reference database, which includes Accumulibacter. Thus, we can just extract
all those reads that map to any of those entries in the reference database
that map to Accumulibacter.

First, we need to extract the reads from the BAM files. This is easily done
in two steps using `samtools`. First, we filter the BAM file for only those
reads that match to a reference of interest. Because we are also dealing with
incomplete genomes that are made up of many contigs, we store all the contig
identifiers in a file `refs.txt`. Then, we can make use of `samtools view`
with a reference filter. The standard command would be `samtools view <refs>
<bamfile>`, where `<refs>` can be multiple references. To read these from the
file, we use `xargs`:

```bash
cat refs.txt | xargs samtools view -b input.bam > filtered.bam
```

Then, we can easily convert these filtered BAM files to FASTQ files:.

```bash
# extract reads as FASTQ
samtools bam2fq -0 pe-null.fastq -1 pe-fwd.fastq -2 pe-rev.fastq \
                -s se.fastq filtered.bam

# combine all reads into a single file and zip to save space
cat pe-null.fastq pe-fwd.fastq pe-rev.fastq se.fastq | gzip > all.fastq.gz

# clean up intermediates to save space
rm pe-null.fastq pe-fwd.fastq pe-rev.fastq se.fastq
```

### 6.2  Make databases for representative genomes

_Bowtie2 databases for these references can be [downloaded
directly](https://www.dropbox.com/s/jkqoub4h8d1puoz/bt2.tar.gz?dl=0)._

The goal here is to identify genomic differences between the Accumulibacter
populations in each granule. We do this by aligning all reads that map to
any Accumulibacter genome to each of the references. This in essence creates
an assembly of the Accumulibacter genomes for each granule, using the known
reference(s) as guides.

Because aligners such as Bowtie2 and BWA only report the best hit, rather
than all hits (unlike BLAST), we cannot directly use the mapping output from
the combined databases, but will re-align all short reads to the reference
genomes.

In our Accumulibacter analysis there are ten reference genomes (click to
download):

- UW1 clade:
  [UW1](/data/Genomes/Accumulibacter/UW1.fasta.gz)
- UW2 clade:
  [UW2](/data/Genomes/Accumulibacter/UW2.fasta.gz) |
  [BA92](/data/Genomes/Accumulibacter/BA92.fasta.gz) |
  [BA93](/data/Genomes/Accumulibacter/BA93.fasta.gz)
- BA91 clade:
  [BA91](/data/Genomes/Accumulibacter/BA91.fasta.gz) |
  [SK01](/data/Genomes/Accumulibacter/SK01.fasta.gz) |
  [SK02](/data/Genomes/Accumulibacter/SK02.fasta.gz)
- BA94 clade:
  [BA94](/data/Genomes/Accumulibacter/BA94.fasta.gz) |
  [SK11](/data/Genomes/Accumulibacter/SK11.fasta.gz) |
  [SK12](/data/Genomes/Accumulibacter/SK12.fasta.gz)

However, in this particular dataset, we only find reads that are assigned to
the clades UW1, UW2, and BA91. Hence, we will only focus on the analysis based
on these three clades.

For the three genomes [UW1](/data/Genomes/Accumulibacter/UW1.fasta.gz),
[UW2](/data/Genomes/Accumulibacter/UW2.fasta.gz),
[BA91](/data/Genomes/Accumulibacter/BA91.fasta.gz) we can build Bowtie2 databases as follows:

```bash
bowtie2-build UW1.fasta.gz UW1
```

This will produce six files `UW1.*.bt2`. We then repeat the same thing for UW2
and BA91.

### 6.3  Map reads to the references

With the extracted reads and mapping databases in hand, we can now re-map all
the reads to each of the genomes. For a single sample and reference genome
(e.g. UW1), the commands for Bowtie2 are:

```bash
# map reads
bowtie2 -x UW1 --no-unal --score-min L,0,-10 \
        -U sample.fastq.gz -S sample-mapped.sam

# sort SAM file and covert to BAM file (BAM is basically compressed and better SAM)
samtools sort sample-mapped.sam > sample-mapped.bam

# create index of BAM file
samtools index sample-mapped.bam
```

We use the following Bowtie2 options:

- `-x`: database
- `--no-unal`: do not keep unmapped reads in the output
- `--score-min L,0,-10`: this sets the scoring threshold. These settings _de
  facto_ turns off thresholding (read more in the Bowtie2 manual)
- `-U`: input reads
- `-S`: output SAM file

Generally, it is never a bad idea to have a quick look at the SAM/BAM files, just to get a feel for what is going on. In principle, SAM files are just text format, so you can just 'look' at them. But for BAM files, you need `samtools view` to translate the output. Also, even if you have SAM files, you might want to use `samtools view` so that you don't see the header information of the file.

```bash
samtools view sample-mapped.bam | less -S
```

We further pipe the output to `less -S` so that it is paged (we can scroll through the output), and also so lines aren't wrapped (`-S`). Each line in a SAM file is a read that is mapped (or potentially unmapped) to a sequence in the database. For example,

```
M01072:70:000000000-ABY0Y:1:2113:11871:6432     16      NC_013194.1     115     6       4M2I16M8I5M4I2M1I2M9I6M1I5M9I9M10I20M1I9M2I10M10I118M1D4M1I173M1I14M4I8M        *       0       0
M01072:70:000000000-ABY0Y:1:2118:16790:18149    16      NC_013194.1     115     6       4M2I16M8I5M4I2M1I2M9I6M1I5M9I9M10I20M1I9M2I9M6I8M1I3M3I108M1D8M1I19M    *       0       0
M01072:70:000000000-ABY0Y:1:1106:17993:5564     16      NC_013194.1     226     42      5M5I4M2I3M3I9M3I74M1D4M1I173M1I7M       *       0       0                                            
M01072:70:000000000-ABY0Y:1:1119:7687:20386     16      NC_013194.1     261     42      60M1D4M1I173M1I14M4I7M2I11M4I9M1I19M6I9M6I8M3I4M4I7M8I13M9I15M  *       0       0                    
M01072:70:000000000-ABY0Y:1:1117:10086:21768    16      NC_013194.1     562     42      10M12I15M9I7M9I5M7I7M1I9M26I14M4I200M1I8M       *       0       0                                    
M01072:70:000000000-ABY0Y:1:1119:11507:7961     16      NC_013194.1     572     7       15M9I7M9I5M7I7M1I9M26I14M4I425M *       0       0                                                    
M01072:70:000000000-ABY0Y:1:2112:12529:8840     16      NC_013194.1     601     42      3M3I21M1I4M1D350M       *       0       0                                                            
M01072:70:000000000-ABY0Y:1:2107:7929:6434      16      NC_013194.1     602     42      4M1I3M1I8M15I6M1I8M11I377M      *       0       0                                                    
M01072:70:000000000-ABY0Y:1:1116:16578:13019    16      NC_013194.1     659     42      4M4I4M5I283M    *       0       0                                                                    
M01072:70:000000000-ABY0Y:1:2118:4369:20230     0       NC_013194.1     699     6       17M2I8M2I4M2I5M4I6M5I2M3I5M3I4M5I10M3I14M2I8M25I43M     *       0       0     
```

Without wanting to understand too much of what is going on here, know only
these facts: The first column is the identifier of the read, and the third
column is the reference it mapped to. The main message is: if there is an
entry, then the read mapped somewhere to the reference.


### 6.3  Aggregate mapped reads (i.e. pileup)

In order to perform comparative genomics, we need to first transform the raw
mapping into a more useful format. One of the basic questions will be: "How to
SNPs (or SNVs) differ between granules?". Thus, what we need is to know the
distribution of base variants for each position of the genome. This is called
a "pile-up" (you can imagine piling up reads on top of the reference). Pileups can be computed using standard tools (such as `samtools`).

```bash
# compute pileup
samtools mpileup -Q 0 -Agf UW1.fasta.gz -t AD sample-mapped.bam > sample-pileup.bcf

# create index
bcftools index ${code}.bcf
```

Let's quickly inspect the output. The top is just header information, but each non-commented line that follows shows the detected bases at each genome position.

```
❯ bcftools view sample-pileup.bcf
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##samtoolsVersion=1.7+htslib-1.7
...

#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  sample-pileup.bcf
NC_013194.1     1       .       A       C,<*>   0       .       DP=1;I16=0,0,1,0,0,0,0,0,0,0,6,36,0,0,0,0;QS=0,1,0;SGB=-0.379885;MQ0F=0 PL:AD   4,3,0,4,3,4:0,1,0
NC_013194.1     2       .       T       C,<*>   0       .       DP=1;I16=0,0,1,0,0,0,0,0,0,0,6,36,0,0,1,1;QS=0,1,0;SGB=-0.379885;MQ0F=0 PL:AD   4,3,0,4,3,4:0,1,0
NC_013194.1     3       .       T       C,<*>   0       .       DP=1;I16=0,0,1,0,0,0,0,0,0,0,6,36,0,0,2,4;QS=0,1,0;SGB=-0.379885;MQ0F=0 PL:AD   4,3,0,4,3,4:0,1,0
...
NC_013194.1     134     .       A       <*>     0       .       DP=6;I16=6,0,0,0,224,8368,0,0,227,9109,0,0,73,1401,0,0;QS=1,0;MQ0F=0    PL:AD   0,18,137:6,0
NC_013194.1     135     .       G       <*>     0       .       DP=6;I16=6,0,0,0,220,8096,0,0,227,9109,0,0,77,1451,0,0;QS=1,0;MQ0F=0    PL:AD   0,18,137:6,0
NC_013194.1     136     .       T       C,<*>   0       .       DP=6;I16=0,0,6,0,0,0,220,8096,0,0,227,9109,0,0,81,1509;QS=0,1,0;VDB=0.0178597;SGB=-0.616816;MQ0F=0      PL:AD   137,18,0,137,18,137:0,6,0
...
```

If we take for example the last line above, the columns are:

- **CHROM**: NC_013194.1, contig name
- **POS**: 4, base position
- ID: .
- **REF**: T, nucleotide in the reference
- **ALT**: G,\* alternate nucleotides detected (a star is a
  placeholder for empty)
- QUAL: 0
- FILTER: .
- **INFO**: DP=1;I16..., additional information

Here, I've highlighted the interesting columns with respect to this analysis
in bold. What this tells us, is that at position 4 of the reference genome,
our sample has alternative nucleotides (in this case G). However, the total
coverage depth (the number of reads that span this genome position is only 1
(DP=1), so it is unclear whether this alternative base is just an error or is
real.

QUAL and FILTER are fields that are used for SNP calling. At this point it is
important to highlight that common SNP callers do not really work as we would
want them to. For one, they assume that the reference is a very good prior for
the *truth*. This need not be the case in our samples, since the reference was
reconstructed from a different data set. This is the main reason why we focus
directly on the SNP profiles.

Finally, the final two columns specify additional information. Here, we have
the fields PL and AD (separated by a :) and then followed by the respective
values. The last line, for example has PL = (137,18,0,137,18,137) and **AD =
(0,6,0)**. This last field is the important one, as is tells how frequent the
different alleles are. For position 136, in this sample we find the nucleotide
`C`, rather than the `T` from the reference.

From this information we can now reconstruct the different distributions of
allele frequencies for each sample.

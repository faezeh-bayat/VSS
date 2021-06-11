# VSS: Variance-stabilized units for sequencing-based genomic signals

## Variance-stabilized signals (VSS) is a signal transformation approach used for eliminating the dependency of data variance from its mean. We generate VSS for sequencing-based genomic signals by learning the empirical relationship between the mean and variance of a given signal data set and producing transformed signals that normalize for this dependence.
## Two main steps in VSS pipeline:
### 1: Training model: Identifying the mean-variance relationship
There are two different options for this step. It either uses the user provided replicates to identify the mean-variance relationship or uses the default trained model. In the latter case, user just needs to provide the untransformed signals.
### 2: Transforming signals: Calculating variance-stabilized signals
Having learned the mean-variance relationship, VSS can be generated using the variance-stabilizing transformation. 



<img src="https://github.com/faezeh-bayat/Variance-stabilized-units-for-sequencing-based-genomic-signals/blob/master/bin/VSS_general_schematic/VSS_schematic.png" width="800"/>

#####################################################################################

## Prerequisites
```
R 4.0.2 (https://www.r-project.org/)
R packages needed
install.packages('bigmemory', 'data.table', 'argparse', 'pracma')

conda install -c macs2
conda install -c bioconda ucsc-bedclip
conda install -c bioconda ucsc-bedgraphtobigwig
conda install -c bioconda ucsc-bigwigtobedgraph

bedtools
######(https://bedtools.readthedocs.io/en/latest/content/installation.html)
```

## Installing VSS pipeline
```
git clone https://github.com/faezeh-bayat/Variance-stabilized-units-for-sequencing-based-genomic-signals.git
```

## How to run VSS pipeline
#### 1. Train the model
##### 1.1 Replicates in bed, bedGraph, bigWig format
```
Rscript VSS.R train rep1 <bed, bedGraph, bigWig> rep2 <bed, bedGraph, bigWig> triandir
               
```

##### 1.2 Tag alignment bam file (Raw read)
```
Rscript VSS.R train_tag tag_alignment_rep1 <bam> tag_alignment_rep2 <bam> --signal "raw" --outdir triandir
              
         
```
##### 1.3 Tag alignment bam file (Fold enrichment and Pvalue)
```
Rscript VSS.R train_tag \
              tag_alignment_rep1 <bam> \
              --fraglen1 <Fragment length for replicate 1> < No need if you choose "raw" signals as output > \
              tag_alignment_rep2 <bam> \
              --fraglen2 <Fragment length for replicate 2> < No need if you choose "raw" signals as output > \
              --chrsz <2-col chromosome sizes file> < No need if you choose "raw" signals as output > \
              --gensz <hs, mm> < No need if you choose "raw" signals as output > \
              --signal <fc, pval, raw, all> 
              --outdir triandir
         
```
#### 2. Transform the signals
##### 2.1 Replicates in bed, bedGraph, bigWig format
```
Rscript VSS.R transform rep1 <bed, bedGraph, bigWig> traindir tranformdir
               
```
##### 2.2 Tag alignment bam file (Raw read)
```
Rscript VSS.R transform_tag rep1 < bam > --signal "raw" traindir tranformdir
               
```
##### 2.3 Tag alignment bam file (Fold enrichment and Pvalue)
```
Rscript VSS.R transform_tag rep1 <bam> \
                            tag_alignment_rep <bam> \
                            --fraglen1 <Fragment length for replicate to be variance stabilized > < No need if you choose "raw" signals as output > \
                            --chrsz <2-col chromosome sizes file> < No need if you choose "raw" signals as output > \
                            --gensz <hs, mm> < No need if you choose "raw" signals as output > \
                            --signal <fc, pval, raw, all> 
                            --model_dir triandir
                            --transdir tranformdir
               
```



#### 3. Variance-stabilized signals will be saved in the tranformdir folder as "Variance_stabilized_signals.bed".


##################################################################################################

#### You can accesee the "Variance-stabilized units for sequencing-based genomic signals" manuscript in:
https://www.biorxiv.org/content/10.1101/2020.01.31.929174v2


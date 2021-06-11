# VSS: Variance-stabilized signals for sequencing-based genomic signals

Variance-stabilized signals (VSS) is a signal transformation approach used for eliminating the dependency of data variance from its mean. We generate VSS for sequencing-based genomic signals by learning the empirical relationship between the mean and variance of a given signal data set and producing transformed signals that normalize for this dependence.
## Two main steps in VSS pipeline:
### 1: Training a model: Identifying the mean-variance relationship
VSS identifies the mean-varaince relationship by using the user provided replicates for an experiment. It needs two replicates for identifying the mentioned relationship.
### 2: Transforming signals: Calculating variance-stabilized signals
Having learned the mean-variance relationship, VSS can be generated using the variance-stabilizing transformation. 



<img src="https://github.com/faezeh-bayat/Variance-stabilized-units-for-sequencing-based-genomic-signals/blob/master/bin/VSS_general_schematic/VSS_schematic.png" width="800"/>



## Prerequisites
```
R 4.0.2 (https://www.r-project.org/)
R packages needed
install.packages('bigmemory', 'data.table', 'argparse', 'pracma')

conda install -c macs2
conda install -c bioconda ucsc-bedclip
conda install -c bioconda ucsc-bedgraphtobigwig
conda install -c bioconda ucsc-bigwigtobedgraph

(You can also directly use the ucsc applications by getting them from "http://hgdownload.soe.ucsc.edu/admin/exe/")
chmod +x ./filePath/utility_name
./filePath/utility_name

bedtools
######(https://bedtools.readthedocs.io/en/latest/content/installation.html)
```

## Installing VSS pipeline
```
git clone https://github.com/faezeh-bayat/VSS.git
```

## How to run VSS pipeline
As mentioned before, VSS has two main steps: Traing a model and transforming the signals. Before traning and transforming signals, VSS pipeline loads the input replicates to convert all data formats to the compatible version used in the pipeline. 

#### 1. Load the data
Input replicates can be in any of bed, bedGraph, bigWig or bam format. In the case that data are in the bam format (tag alignment data), you have multiple options. You can either convert the bam file to raw signals or you can convert them to any of "Fold enrichment (fc)" or "p-value (pval)" signals. We seperate these two conditions as you need to provide more arguments to pipeline to convert the bam file to either of fc or pval signals. We use ENCODE's default parameters for calculating the fc/pval signals.

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


# TensorQTL workshop for the mQTL mapping using covariates

Before starting with the workshop, you should have [PLINK](https://www.cog-genomics.org/plink/1.9/), python3, [TensorQTL and dependencies](https://github.com/broadinstitute/tensorqtl).

The inputs for TensorQTL in the case of performing cis-mQTLs with covariates are: 
- PLINK binary file (.bim, .bed, .fam) format for the genotype. 
- BED UCSC file for the phenotype 
- Text file for the covariates. 

## Access to virtual machine

For accessing the virtual machine that you installed, you need to type the following command on your command shell:

`ssh user@127.0.0.1 -p 5679`

## Genotype input

The input for TensorQTL of the genotype data must be a PLINK binary file format, which is a set of three different files: [.bim](https://www.cog-genomics.org/plink/1.9/formats#bim), containing the information of the SNPs, [.fam](https://www.cog-genomics.org/plink/1.9/formats#fam), with the information of the samples, and [.bed](https://www.cog-genomics.org/plink/1.9/formats#bed), a binary file with the genotype of the samples. But in some cases, the genotype raw data comes within a [Variant Call Format file](https://samtools.github.io/hts-specs/VCFv4.2.pdf). Therefore, with the PLINK tools, we will need to transform it: 

`plink --vcf {input name} --make-bed --out {output name}`

If you wanna have a look at the files: 
```
#to see info of the samples
head chr22.fam 
#to see info fo the SNPs
head chr22.bim 
#to see the quantitative genotype
hexdump -x chr22.bed | head 
```  

## Phenotype input

In this case, for the phenotype, we need a [BED UCSC file](https://genome.ucsc.edu/FAQ/FAQformat.html#format1). Coming from the Quality Control performed into the data using minfi Bioconductor's package or other tools, we will get an [ExpressionSet](https://www.bioconductor.org/packages/release/bioc/vignettes/Biobase/inst/doc/ExpressionSetIntroduction.pdf) which contains the phenoData, the annotation of the probes and the beta values. By R we will create the template of the BED file by writing a text file with the chromosome, the start, the end, the ID and the beta values of the probes per sample. Once got it, with the following command we will sort it by the genomic coordinates: 

`(head -n1 chr22_bed_file.txt && sort -k1,1V -k2,2n -k3,3n <(tail -n+2 chr22_bed_file.txt)) > chr22_phenotype.bed`

Once we got the file sorted, we will zip it and index it:

```
bgzip chr22_phenotype.bed

tabix -p bed chr22_phenotype.bed.gz
```

At the end the file will look like: 

```
#Chr start end ID sample1 sample2 sample3 sample4 
chr1 173863 173864 ENSG123 -0.50 0.82 -0.71 0.83
chr1 685395 685396 ENSG456 -1.13 1.18 -0.03 0.11
chr1 700304 700305 ENSG789 -1.18 1.32 -0.36 1.26
```
**Notice that this example doesn't belong to methylation data, but to expression data, as we already know the beta values are between 0 and 1,** to see it: 

`zcat chr22_phenotype.bed.gz`

## Covariates input

For the covariates, we will use a text file, in which the first line corresponds to the sample names and the following ones to the covariates we wish, e.g. sex, Principal Components, batch... At the end, the file should look like this: 

```
id sample1 sample2 sample3 sample4
PC1 -0.02 0.14 0.16 -0.02
PC2 0.01 0.11 0.10 0.01
PC3 0.03 0.05 0.08 0.07
```
If you do `head chr22_covariables.txt` you will see that the covariates used are the sex of the samples and the Prinicipal Components from 1 to 5. 

## Run TensorQTL for cis-mQTL mapping with covariates

The first step is to open python3 by typing it on the command line: 

`python3` 

The next step is to import the necessary libraries: 

```
import pandas as pd
import torch
import tensorqtl
from tensorqtl import genotypeio, cis, trans
print('PyTorch {}'.format(torch.__version__))
print('Pandas {}'.format(pd.__version__))
```

To read the files we will save the path to the following covariates: 

```
plink_prefix_path = '/home/user/workshop/chr22'
expression_bed = '/home/user/workshop/chr22_phenotype.bed.gz'
covariates_file = '/home/user/workshop/chr22_covariables.txt'
```

Upload the files: 

```
# load phenotypes and covariates
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T

# PLINK reader for genotypes
pr = genotypeio.PlinkReader(plink_prefix_path)
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]
```

We need to have the same order of the samples between the covariates and the phenotype: 

`phenotype_df = phenotype_df.reindex(sorted(phenotype_df.columns), axis=1)`

Run TensorQTL: 

`cis_df = cis.map_cis(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df=covariates_df, seed=123456789)`

- genotype_df contains the quantitative genotype.
- variant_df contains the information regarding the SNPs. 
- phenotype_df contains the beta values. 
- phenotype_pos_df contains the annotation of the probes. 
- covariates_df contains the covariates choosen. 
- seed, we use it to be able to replicate the experiment. 

To save the results in a text file and be able to analyse them on another platform: 

`cis_df.to_csv('/home/user/workshop/cis_tensorQTL_chr22.txt', header=True, index=True, sep='\t')`

In case that you haven't been able to run TensorQTL, you will find available the results in this [link](https://ehubox.ehu.eus/s/WNFxR97nzsNELoo). 

Finally, to correct for multiple testing, we will use R: 

```
tensor <- read.table("/home/user/workshop/cis_tensorQTL_chr22.txt",sep ="\t", header = T)
tensor$bonferroni = p.adjust(tensor$pval_beta, method="bonferroni")
```

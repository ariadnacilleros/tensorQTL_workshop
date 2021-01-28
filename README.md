# TensorQTL workshop for the mQTL mapping using covariates

Before starting with the workshop, you should have [PLINK](https://www.cog-genomics.org/plink/1.9/), python3, [TensorQTL and dependencies](https://github.com/broadinstitute/tensorqtl).

The inputs for TensorQTL in the case of performing cis-mQTLs with covariates are: 
- PLINK binary file (.bim, .bed, .fam) format for the genotype. 
- BED UCSC file for the phenotype 
- Text file for the covariates. 

## Genotype input

The input for TensorQTL of the genotype data must be a PLINK binary file format, which is a set of three different files: [.bim](https://www.cog-genomics.org/plink/1.9/formats#bim), containing the information of the SNPs, [.fam](https://www.cog-genomics.org/plink/1.9/formats#fam), with the information of the samples, and [.bed](https://www.cog-genomics.org/plink/1.9/formats#bed), a binary file with the genotype of the samples. But on some cases, the genotype raw data comes within a [Variant Call Format file](https://samtools.github.io/hts-specs/VCFv4.2.pdf). Therefore, with the PLINK tools, we will need to transform it: 

`plink --vcf {input name} --make-bed --out {output name}`

## Phenotype input

On this case, for the phenotype we need a [BED UCSC file](https://genome.ucsc.edu/FAQ/FAQformat.html#format1). Coming from the Quality Control performed into the data using minfi Bioconductor's package or another tools, we will get an [ExpressionSet](https://www.bioconductor.org/packages/release/bioc/vignettes/Biobase/inst/doc/ExpressionSetIntroduction.pdf) which contains the phenoData, the annotation of the porbes and the beta values. By R we will create the template of the BED file by writing a text file with the chromosome, the start, the end, the ID and the beta values of the probes per sample. Once got it, with the following command we will sort it by the genomic coordinates: 

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
**Notice that this example doesn't belong to methylation data, but to expression data, as we already know the beta values are between 0 and 1.** 

## Covariates input

For the covariates we will use a text file, in which the first line corresponds to the sample names and the following ones to the covariates we wish, e.g. sex, Principal Components, batch... At the end the file should look like this: 

```
id sample1 sample2 sample3 sample4
PC1 -0.02 0.14 0.16 -0.02
PC2 0.01 0.11 0.10 0.01
PC3 0.03 0.05 0.08 0.07
```

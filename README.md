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

## Covariates input


# BIOE 582 Project 1

Coursework in computational genomics completed as part of my Master's of Engineering in Bioinformatics at the University of Illinois Chicago.

## Part 1: Affymetrix Array Preprocessing

The assignment introduces preprocessing Affymetrix oligonucleotide array probe-level data with Bioconductor's `affy` package. Its dilution dataset includes human liver and CNS-derived RNA hybridized to HGU95A arrays at known proportions and concentrations, with replicate samples processed on different scanners. These known inputs make it possible to examine whether measured expression tracks RNA abundance and whether replicates agree after preprocessing.

[Read the Part 1 assignment (PDF)](bioe582pj1b_mcgreevy.pdf)

## Part 2: Microarray Preprocessing and Annotation

The R script accepts a sample-and-phenotype metadata file and an output directory, then processes four Affymetrix CEL files split between estrogen receptor-positive (ER+) and estrogen receptor-negative (ER-) samples. It applies RMA background correction and normalization, and saves side-by-side boxplots comparing raw and normalized measurements. MAS 5.0 detection calls are used to filter probe sets that are not present in at least one sample of each phenotype; Affymetrix control probes are also removed. The script reports filtered probes and writes the remaining probes' average expression by phenotype, annotated with Entrez IDs and gene symbols.

[View the Part 2 R script](bioe582pj1b_mcgreevy_main.R)


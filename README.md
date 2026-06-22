---
editor_options: 
  markdown: 
    wrap: 72
---

# Do body size mean and variance share the same genetic architecture?

Internship Edouard Courret 2026

This is all the research for the subject: "Do body size mean and
variance share the same genetic architecture?". The aim is to find
common QTLs associated with mean and variance in *C.elegans*.

For that, we studied the mean body size and intra-genotypic variance
with two metrics: variance corrected by the mean and coefficient of
variation.

Here you can find the data, plot and script for analysis.

### Summary of the script:

 - 1: Clean the data before analysis

\- 2: BLUEs of mean

\- 3,3bis: BLUEs of variance and coefficient of variation (CV)

\- 4: Preparing Files for GEMMA

\- 5,6,7: Analysis of GEMMA results for mean, variance and CV

\- 8: Preparing Files for GridLMM

\- 9: Preparing Files for permGWAS

\- 10: function to run gridLMM (NGM)

\- 11: Preparation of GridLMM results and pruning

\- 12: function to run gridLMM (NaCl)

\- 13: GWA analysis of gridLMM results

\- 14: Estimations of BLUEs in NGM and NaCl, correlations mean-variance,
distributions of phenotypic values

\- 15: Estimations of heritability

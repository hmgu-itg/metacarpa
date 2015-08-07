#Metacarpa: META-analysis in C++ Accounting for Relatedness, using arbitrary Precision Arithmetic

## Description
This is the repository where the code for METACARPA lives. 

## Background
As open data and data sharing policies increase among the scientific communities, more and more studies are expected to share their results in the coming years. This opens exciting perspectives for meta-analysis studies, which can aggregate such results in order to boost power and discover potential new genetic associations. When performing meta-analysis, particular care needs to be given to sample selection, because meta-analysed studies are supposed to contain independent samples only.However, due to privacy policies, researchers often do not publish raw genotype data, or if they do, they scramble sample IDs so that no genotype can be traced back to the actual person.
METACARPA is designed for meta-analysing genetic association studies with overlapping or related samples, when details of the overlap or relatedness are unknown. 

It implements and expands a method first described by [Province and Borecki](www.ncbi.nlm.nih.gov/pmc/articles/PMC3773990/).
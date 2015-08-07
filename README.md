#Metacarpa: META-analysis in C++ Accounting for Relatedness, using arbitrary Precision Arithmetic

## Description
This is the repository where the code for METACARPA lives. 

## Background
As open data and data sharing policies increase among the scientific communities, more and more studies are expected to share their results in the coming years. This opens exciting perspectives for meta-analysis studies, which can aggregate such results in order to boost power and discover potential new genetic associations. When performing meta-analysis, particular care needs to be given to sample selection, because meta-analysed studies are supposed to contain independent samples only.However, due to privacy policies, researchers often do not publish raw genotype data, or if they do, they scramble sample IDs so that no genotype can be traced back to the actual person.

METACARPA is designed for meta-analysing genetic association studies with overlapping or related samples, when details of the overlap or relatedness are unknown. It implements and expands a method first described by [Province and Borecki](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3773990/). In this article, the authors describe how to meta-analyse p-values, however what researchers are most interested in are effect sizes, i.e. the combined effect of a particular SNP across all studies. [Lin and Sullivan](http://www.ncbi.nlm.nih.gov/pubmed/20004761) describe how this can be done when the degree of overlap is known; METACARPA combines both these methods.

## Installation

### Binaries

METACARPA is compiled as a static Linux x64 executable. Running it should be as simple as downloading [this executable](linux_64_static_bin/metacapa2) and running it. 

> The program is still under development and deployment has not yet been tested. Please email the author if you encounter any problems.

### Compilation

The src directory contains a single source file and Makefile. The [Makefile](src/Makefile) is written for a specific environment and is not deployment-ready. If you plan to try and compile METACARPA, you will need to change at least the first three lines of the Makefile to add custom compiler paths. Please note that you will also need a copy of the [Boost](http://www.boost.org) mathematical libraries on your machine (the binaries were compiled using Boost 1.55.0).

## Running

Here is a short summary of METACARPA's options. 
```
METACARPA ( μ - 🐟  ): Meta-analysis in C++ Accounting for Relatedness using arbitrary Precision Arithmetic.
===========================================================================================================

        NB: All arguments mandatory except column arguments.
        MATACARPA currently supports only one header line.

Options description :
  --help                This help message.
  -I [ --input ] arg    Input file.
  -O [ --output ] arg   Output file.
  -t [ --sep ] arg      Input field separator. Don't forget to quote if 
                        necessary. Output field separator is always \t.
  -c [ --chr-col ] arg  1-based p-value column number.
  -q [ --pos-col ] arg  1-based position column number.
  -a [ --all-col ] arg  1-based column number for effect or reference allele.
  -r [ --rsid-col ] arg 1-based column number for RSID or any other column that
                        you want to keep.
  -p [ --pval-col ] arg 1-based p-value column number.
  -b [ --beta-col ] arg 1-based beta column number.
  -s [ --se-col ] arg   1-based beta-SE column number.
  -n [ --size-col ] arg 1-based sample size column (if absent, sample sizes 
                        will be assumed constant and should be appended to 
                        input file names using a comma : -I 
                        [FILENAME],[SAMPLE_SIZE]).
  -i [ --id-col ] arg   1-based ID column number (must be unique - e.g. 
                        chr:pos-A1-A2). If absent, chr:pos will be used.
  -m [ --matrix ] arg   Path to a METACARPA-generated correlation matrix array.

```
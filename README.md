#Metacarpa: META-analysis in C++ Accounting for Relatedness, using arbitrary Precision Arithmetic

## Description
This is the repository where the code for METACARPA lives. 

## Background
As open data and data sharing policies increase among the scientific communities, more and more studies are expected to share their results in the coming years. This opens exciting perspectives for meta-analysis studies, which can aggregate such results in order to boost power and discover potential new genetic associations. When performing meta-analysis, particular care needs to be given to sample selection, because meta-analysed studies are supposed to contain independent samples only.However, due to privacy policies, researchers often do not publish raw genotype data, or if they do, they scramble sample IDs so that no genotype can be traced back to the actual person.

METACARPA is designed for meta-analysing genetic association studies with overlapping or related samples, when details of the overlap or relatedness are unknown. It implements and expands a method first described by [Province and Borecki](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3773990/). In this article, the authors describe how to meta-analyse p-values, however what researchers are most interested in are effect sizes, i.e. the combined effect of a particular SNP across all studies. [Lin and Sullivan](http://www.ncbi.nlm.nih.gov/pubmed/20004761) describe how this can be done when the degree of overlap is known; METACARPA combines both these methods.

## Installation

### Binaries

METACARPA is compiled as a static Linux x64 executable. Running it should be as simple as downloading [this executable](https://bitbucket.org/agilly/metacarpa/downloads/metacarpa) and running it. 

> The program is still under development and deployment has not yet been tested extensively. Please email the author if you encounter any problems.

### Compilation

The src directory contains a single source file and Makefile. The [Makefile](src/Makefile) is written for a compiler that supports C++ 11. Please note that you will also need a copy of the [Boost](http://www.boost.org) mathematical libraries on your machine (the binaries were compiled using Boost 1.55.0).

## Program Options

Here is a short summary of METACARPA's options. 
```
METACARPA ( μ - 🐟  ): Meta-analysis in C++ Accounting for Relatedness using arbitrary Precision Arithmetic.
===========================================================================================================

        Usage : /nfs/users/nfs_a/ag15/metacarpa-git/src/metacarpa -I infile1[,size] [-I infile2[,size] ...] -O outfile --chr-col int --pos_col int --a1-col int--a2-col int --pval-col int --beta-col int --se-col int --af-col int [--size-col int] [--sep char] [--id-col int] [-m matrix_file] [-x] [-d]

NB:     METACARPA currently supports only one header line in input files, which is ignored.

Options description :
  --help                This help message.
  -I [ --input ] arg    (mandatory) Input file.
  -O [ --output ] arg   (mandatory) Output file.
  -c [ --chr-col ] arg  (mandatory) 1-based chromosome column number.
  -q [ --pos-col ] arg  (mandatory) 1-based position column number.
  -u [ --a1-col ] arg   (mandatory) 1-based column number for effect or
                        reference allele.
  -v [ --a2-col ] arg   (mandatory) 1-based column number for other allele.
  -p [ --pval-col ] arg (mandatory) 1-based p-value column number.
  -b [ --beta-col ] arg (mandatory) 1-based beta column number.
  -s [ --se-col ] arg   (mandatory) 1-based beta-SE column number.
  -a [ --af-col ] arg   (mandatory) 1-based effect allele frequency.
  -n [ --size-col ] arg 1-based sample size column (if absent, sample sizes
                        will be assumed constant and should be appended to
                        input file names using a comma : -I infile,size).
  -t [ --sep ] arg      Input field separator. Don't forget to quote if
                        necessary. Output field separator is always \t.
  -i [ --id-col ] arg   1-based RSID column number (should be unique - e.g.
                        chr:pos-A1-A2). If absent, chr:pos will be used.
  -m [ --matrix ] arg   Path to a METACARPA-generated correlation matrix array.
  -x [ --stop ]         Stop METACARPA after generating the matrix.
  -d [ --debug ]        Toggles an extremely verbose output, for debugging
                        purposes only.
```

## Data preparation

Although METACARPA scales well to large numbers of variants, you can LD-prune or randomly thin down your SNPs to accelerate the calculation of study correlations. This can be done, for example, using [Plink](http://cog-genomics.org/plink2/). For more detailed estimations of minimum number of SNPs needed to get correct results, please refer to the [Metacarpa simulation report](bitbucket.org/agilly/metacarpa-simulation).

## File formats

METACARPA supports any fixed-delimiter input file that has the column it requires. You should specify the column indices that correspond to each of these columns using the command line arguments. For example:

`./metacapa2 -I file1.assoc,1500 -I file2.assoc,2300 -O out.meta_analysis.txt -t '\t' -c 1 -q 2 -r 3 -i 3 -p 4 -b 5 -s 6`

will run the meta-analysis on two studies with $1500$ and $2300$ samples respectively. If your input format contains the number of individuals genotyped at a particular position, remove the commas and sample sizes after the file name and replace by `-n sample_col`, where `sample_col` is the column number that contains this information.

> **WARNING: ** The input files should be sorted according to chromosome and position.

## Field separators

On some systems, it can be complicated to pass escaped characters such as tabs (`\t`) as parameters on the command line. In Bash, for example, you can circumvent this by telling the shell not to interpret the character, prefixing it with a dollar sign, like so: `$'\t'`.

## Output format
METACARPA produces an output with the following columns:
* **rsid** : RSid, this should be equal to the ID column (argument -i) in the input.
* **chr:pos** : SNPid.
* **effect_allele** : Effect allele.
* **neffect_allele** : Non-effect allele.
* **effect_allele_frequency** : Total effect allele frequency.
* **effects** : Summary of directions of effect across files.
* **beta** : Meta-analysis Effect.
* **se** : Meta-analysis Standard error.
* **z** : Meta-analysis Z-score.
* **zse** : Standard error of the z-score.
* **p_wald** : Corrected p-value from the effect-size based meta-analysis. This is the preferred p-value to use.
* **p_corrected** : Corrected p-value from the p-value based meta-analysis. This is the adapted p-value from Province and Borecki (2013).
* **p_stouffer** : Stouffer p-value. This is an ** * uncorrected * ** p-value and should be used only for comparison purposes.

## Allele formatting
METACARPA expects alleles in capital letters, small letters will be considered as different alleles. METACARPA implements allele flipping, but does not handle strand issues.

## Multiple runs

METACARPA runs in two steps:

* calculation of the correlation matrix
* SNP per SNP meta-analysis.

If you run meta-analyses for different traits, if you chunk your data or generally if you want to meta-analyse the same populations more than once, you can reuse the correlation matrix by using the `-m` argument. By default, METACARPA will always write the matrices it calculates to disk, they should appear as small `*.matrix` files in the running directory. If you point the `-m` argument to one of these files, METACARPA will skip calculation of the matrix and jump directly to the meta-analysis. The `--stop` or `-x` argument tells METACARPA to only calculate the matrix, without performing any single-point analysis. 

> If you decide to chunk your data, please ensure that you calculate the correlation matrix on the genome-wide, un-chunked set. You can then use that matrix on your chunked data for meta-analysis.

## Scalability

The matrix calculation step is the most costly part of the algorithm, but it has been designed to scale well to very large datasets such as sequencing data. A memory assignment of 1Gb should be more than enough for most cases including whole-genome sequence data.

## Performance and accuracy

For detailed analyses of this software's performance and comparisons with other methods, please refer to the [Metacarpa simulation report](bitbucket.org/agilly/metacarpa-simulation).
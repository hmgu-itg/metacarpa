# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

METACARPA (META-analysis in C++ Accounting for Relatedness, using arbitrary Precision Arithmetic) is a C++11 bioinformatics tool for meta-analysing genetic association studies with overlapping or related samples when the degree of overlap is unknown. It combines the p-value correction method of Province and Borecki (2013) with the effect-size method of Lin and Sullivan (2009).

## Build

The project is a single C++ source file compiled via Make. Boost 1.60.0+ is required (program_options, serialization, multiprecision, uBLAS, math).

```bash
# Build (produces static Linux x64 binary at src/metacarpa)
make -C src/

# The Makefile IDIR and LDIR variables must point to your Boost headers and libs
```

Compiler flags: `-O3 -std=c++11 -static`. Links against: `boost_program_options`, `boost_serialization`, `pthread`.

## Architecture

The entire implementation lives in a single file: `src/metacarpa.cpp` (~1850 lines).

**Two-phase algorithm:**
1. **Matrix calculation** — Reads all input files, identifies overlapping variant positions across studies, and computes a variance-covariance correlation matrix using tetrachoric correlation of binarized p-values. The matrix is serialized to disk (`.matrix.txt`) for reuse across runs.
2. **Single-point meta-analysis** — Re-reads input files and performs per-variant meta-analysis using the precomputed correlation matrix, producing corrected effect sizes, Z-scores, and p-values.

**Key components in `metacarpa.cpp`:**
- `studies_correlation` class — Stores and manages correlation matrices indexed by study "mask" (bitmask of which studies contribute to a variant). Handles Boost serialization (XML/text).
- `count_data` struct (nested in `studies_correlation`) — 2x2 contingency tables for tetrachoric correlation.
- `output_record` struct — Holds per-variant meta-analysis results.
- `meta_analyse()` — Core computation: inverts correlation submatrix, computes corrected Z-scores and effect sizes.
- `z_transform()` / `z_transform_fess()` — P-value to Z-score conversion with automatic fallback to `cpp_dec_float<200>` arbitrary precision for extremely small p-values (< ~1e-29).
- `InvertMatrix()` — LU-decomposition matrix inversion via Boost uBLAS.
- `tetrachoric()` — Estimates correlation between studies from binary-transformed p-value vectors.

**Data flow:** Parse CLI args → open sorted input files → first pass builds correlation matrix → serialize matrix → second pass performs per-variant meta-analysis → write tab-delimited output.

## Key Constraints

- Input files must be sorted by chromosome then position.
- Alleles are case-sensitive (expects capitals); allele flipping is supported but strand flipping is not.
- Sample sizes can be per-variant (via `--size-col`) or constant (appended to filename: `-I file,N`).
- The correlation matrix can be precomputed with `-x` (stop after matrix) and reused with `-m` for different traits or chunked data.

## No Test Suite

There is no automated test framework. Validation is done externally via the [metacarpa-simulation](https://bitbucket.org/agilly/metacarpa-simulation) project.

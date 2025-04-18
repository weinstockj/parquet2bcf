# Introduction

parquet2bcf is a utility that converts variant genetic data from Parquet format to BCF (Binary Call Format).

## Overview

This tool reads sparse genotype data from a Parquet file with columns for chromosome, position, reference allele, alternate allele, and sample ID, then writes it to a BCF file that can be used with bioinformatics tools.

## Features

- Fast conversion from Parquet to BCF
- Support for sparse genotype representation
- Multi-threaded BCF writing
- Configurable input/output paths
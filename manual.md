---
layout: page
title: "Manual"
description: ""
group: navigation
---
{% include JB/setup %}

Typing `kallisto` produces a list of usage options, which are:

    Usage: kallisto <CMD> [arguments] ..

    Where <CMD> can be one of:

        index         Builds a kallisto index 
        quant         Runs the quantification algorithm 
        h5dump        Converts HDF5-formatted results to plaintext
        version       Prints version information

    Running kallisto <CMD> without arguments prints usage information for <CMD>

The four usage commands are:

#### index

`kallisto index` builds an index from a Fasta formatted file of target sequences. The arguments for the index command are:

    Usage: kallisto index [arguments] FASTA-file

    Required argument:
    -i, --index=STRING          Filename for the kallisto index to be constructed 
     
    Optional argument:
    -k, --kmer-size=INT         k-mer (odd) length (default: 31, max value: 31)

The Fasta file supplied can be either in plaintext or gzipped format. 

#### quant

`kallisto quant` runs the quantification algorithm. The arguments for the quant command are:

    Usage: kallisto quant [arguments] FASTQ-files
    
    Required arguments:
    -i, --index=INT               Filename for the kallisto index to be used for
                                  quantification
    -o, --output-dir=STRING       Directory to write output to
    
    Optional arguments:
    -l, --fragment-length=DOUBLE  Estimated average fragment length
                                  (default: value is estimated from the input data)
    -b, --bootstrap-samples=INT   Number of bootstrap samples (default: 0)
        --seed=INT                Seed for the bootstrap sampling (default: 42)
        --plaintext               Output plaintext instead of HDF5

#### h5dump

`kallisto h5dump` converts [HDF5](https://www.hdfgroup.org/HDF5/whatishdf5.html)-formatted results to plaintext. The arguments for the h5dump command are:

    Usage:  kallisto h5dump [arguments] abundance.h5

    Required argument:
    -o, --output-dir=STRING       Directory to write output to

#### version

`kallisto version` displays the current version of the software.
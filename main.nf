#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process celluster {

    input:
    // list of cluster files
    // The orignal file that was used for clustering.

    script:
    '''
    python3 celluster.py
    '''
}

workflow {
    celluster
}
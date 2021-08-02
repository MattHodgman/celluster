#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process celluster {

    input:
    path input_files // cluster labels from different methods
    path data // original data file used for clustering

    output:
    path 'exemplar-001-iter-consensus.csv'
    path 'exemplar-001-iter-outliers.csv'

    script:
    """
    python3 celluster.py -i $input_files*.csv -c CellID -d $data -v -n exemplar-001 -r 1
    """
}

workflow {
    input_files = channel.fromPath('data/*cells.csv')
    data = channel.fromPath('data.unmicst-exemplar-001.csv')
    celluster(input_files, data)
}
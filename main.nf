#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process celluster {

    input:
    path input_files // cluster labels from different methods
    path data // original data file used for clustering

    script:
    """
    python3 celluster.py -i $input_files*.csv -c CellID -d $data -v -n test -r 1
    """
}

workflow {
    input_files = channel.fromPath('/cluster_labels/*.csv')
    data = channel.fromPath('data.csv')
    celluster(input_files, data)
}
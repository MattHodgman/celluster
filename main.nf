#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.iter = 1
params.name = 'test'
params.id = 'CellID'


process fastpg {
    publishDir "./data/fastpg/", mode: 'move'

    input:
    path data

    output:
    path '*.csv'

    script:
    """
    python3 /app/cluster.py -i $data -c
    """
}

process scanpy {
    publishDir "./data/scanpy/", mode: 'move'

    input:
    path data
    
    output:
    path '*.csv'
    
    script: 
    """
    python3 /app/cluster.py -i $data -c
    """
}

process flowsom {
    publishDir "./data/flowsom/", mode: 'move'

    input:
    path data

    output:
    path '*.csv'
    
    script:
    """
    python3 /app/cluster.py -i $data -c
    """
}

process celluster {

    publishDir ".", mode: 'move'

    input:
    path input_files // cluster labels from different methods
    path clean_data // original data file used for clustering

    output:
    path '*.csv'

    script:
    """
    python3 /app/celluster.py -i $input_files -c $params.id -d $clean_data -v -r $params.iter
    """
}

workflow {
    // input raw data file
    data = channel.fromPath('data/unmicst-exemplar-001.csv')

    // cluster with different methods
    fastpg(data)
    scanpy(data)
    flowsom(data)

    // get output cell cluster label files
    input_files = channel.fromPath('data/*/*-cells.csv')
    clean_data = channel.fromPath('data/*/*-clean.csv').first()

    // run consensus clustering
    celluster(input_files, clean_data)
}
#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.iter = 1
params.name = 'test'
params.id = 'CellID'

process fastpg {

}

process scanpy {

}

process flowsom {
    
}

process celluster {

    publishDir ".", mode: 'symlink'

    input:
    path input_files // cluster labels from different methods
    path data // original data file used for clustering

    output:
    path '*.csv'

    script:
    """
    celluster.py -i $input_files -c $params.id -d $data -v -n $params.name -r $params.iter
    """
}

workflow {
    input_files = channel.fromPath('data/*cells.csv').toList()
    data = channel.fromPath('data/unmicst-exemplar-001.csv')
    celluster(input_files, data)
}
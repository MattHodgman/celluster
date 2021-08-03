#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.id = 'CellID'
params.input = ''

process fastpg {

    input:
    tuple val(dataName), file(dataFile)

    output:
    path '*-cells.csv', emit: cells

    script:
    """
    python3 /app/cluster.py -i ${dataFile} -c
    """
}

process scanpy {

    input:
    tuple val(dataName), file(dataFile)
    
    output:
    path '*-cells.csv', emit: cells
    
    script: 
    """
    python3 /app/cluster.py -i ${dataFile} -c
    """
}

process flowsom {

    input:
    tuple val(dataName), file(dataFile)

    output:
    path '*-cells.csv', emit: cells
    
    script:
    """
    python3 /app/cluster.py -i ${dataFile} -c
    """
}

process celluster {

    publishDir ".", mode: 'move'

    input:
    file '*.csv' // cluster labels from different methods
    tuple val(dataName), file(dataFile) // original data file used for clustering and name
    val iter

    output:
    path '*-iter-consensus.csv', emit: consensus
    path '*-iter-outliers.csv', emit: outliers

    script:
    """
    python3 /app/celluster.py -i *.csv -c $params.id -d ${dataFile} -v -r ${iter} -n ${dataName}
    """
}

def len(f) {
    f.subscribe { println it.readLines().size() }
}

workflow {
    // input raw data file
    data = channel
                .fromPath(params.input)
                .map { file -> tuple(file.baseName, file) }

    iter = channel.value(1)
    iter = iter + 1

    // cluster with different methods
    fastpg(data)
    scanpy(data)
    flowsom(data)

    // get output cell cluster label files
    input_files = scanpy.out.cells
                                .mix(fastpg.out.cells, flowsom.out.cells)
                                .toList()

    // run consensus clustering
    celluster(input_files, data, iter)
    len(celluster.out.outliers)
}
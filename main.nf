#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.id = 'CellID'
params.input = ''
params.iter = 1
params.clusters = 'NO_FILE'
params.name = ''

process fastpg {

    publishDir 'fastpg'

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

    publishDir 'scanpy'

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

    publishDir 'flowsom'

    input:
    tuple val(dataName), file(dataFile)

    output:
    path '*-cells.csv', emit: cells
    
    script:
    """
    python3 /app/cluster.py -i ${dataFile} -c -n 15
    """
}

process celluster {

    publishDir 'consensus', mode: 'move'

    input:
    file '*.csv' // cluster labels from different methods
    tuple val(dataName), file(dataFile) // original data file used for clustering and name
    path consensus

    output:
    path '*-iter-consensus.csv'
    path '*-iter-outliers.csv'
    path '*-running-consensus.csv' optional true

    script:
    if (params.iter == 1) {
        """
        python3 /app/celluster.py -i *.csv -c $params.id -d ${dataFile} -r $params.iter -n ${dataName} -z
        """
    } else if (params.iter > 1) {
        """
        python3 /app/celluster.py -i *.csv -c $params.id -d ${dataFile} -r $params.iter -n $params.name -k $consensus -z
        """
    }
}


workflow {
    // input raw data file
    data = channel
                .fromPath(params.input)
                .map { file -> tuple(file.baseName, file) }

    // make method output file dir
    file('./fastpg').mkdir()
    file('./scanpy').mkdir()
    file('./flowsom').mkdir()
    file('./consensus').mkdir()

    // cluster with different methods
    fastpg(data)
    scanpy(data)
    flowsom(data)

    // get output cell cluster label files
    input_files = scanpy.out.cells
                                .mix(fastpg.out.cells, flowsom.out.cells)
                                .toList()

    consensus = channel.fromPath(params.clusters)

    // run consensus clustering
    celluster(input_files, data, consensus)
    
    // celluster.out.outliers.subscribe { 
    //     if ( it.readLines().size() > 100) {
    //         println ">100 outliers!"
    //     }
    // }
}
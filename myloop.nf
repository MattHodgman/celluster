#!/usr/bin/env nextflow

nextflow.enable.dsl=2

import nextflow.extension.CH
import static nextflow.extension.DataflowHelper.newOperator
import groovyx.gpars.dataflow.operator.ChainWithClosure
import groovyx.gpars.dataflow.operator.CopyChannelsClosure

params.id = 'CellID'
params.input = ''

def attach(source, target) {
    newOperator([source.createReadChannel()], [target],
                new ChainWithClosure(new CopyChannelsClosure()))
}

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
    python3 /app/celluster.py -i *.csv -c $params.id -d ${dataFile} -n ${dataName} -r ${iter}
    """
}

workflow run_celluster {

  take:
    input_ch
    iter

  main:

    iter.view()

    // input raw data file
    data = input_ch
                .map { file -> tuple(file.baseName, file) }

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

    celluster.out.outliers.view()

  emit:
    celluster.out.outliers
}

iter = channel.value(0)

workflow {
    condition = { it.readLines().size()<100 }
    feedback_ch = CH.create()
    input_ch = Channel.fromPath(params.input, checkIfExists:true)
                      .mix( feedback_ch.until(condition) )

    input_ch.view()
    iter++

    cell_ch = run_celluster(input_ch, iter)

    attach(cell_ch, feedback_ch)
}
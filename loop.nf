#!/usr/bin/env nextflow

nextflow.enable.dsl=2

import nextflow.extension.CH
import static nextflow.extension.DataflowHelper.newOperator
import groovyx.gpars.dataflow.operator.ChainWithClosure
import groovyx.gpars.dataflow.operator.CopyChannelsClosure

def attach(source, target) {
    newOperator([source.createReadChannel()], [target],
                new ChainWithClosure(new CopyChannelsClosure()))
}

process foo {
    input:
      path x
    output:
      path 'foo.txt'
    script:
    """
    cat $x > foo.txt
    """
}

process bar {
    input:
      path x
    output:
      path 'bar.txt', emit: feedback
      path 'bar.txt', emit: result
    script:
    """
    cat $x > bar.txt
    echo World >> bar.txt
    """
}

params.input = "$baseDir/hello.txt"
workflow {
    condition = { it.readLines().size()>4 }
    feedback_ch = CH.create()
    input_ch = Channel.fromPath(params.input, checkIfExists:true)
                      .mix( feedback_ch.until(condition) )
    foo_ch = foo(input_ch)
    bar_ch = bar(foo_ch)
    attach(bar_ch.feedback, feedback_ch)
    bar_ch.result.last().view { "Result:\n${it.text.indent(' ')}" }
}
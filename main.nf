#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { REALIGNMENT } from './workflows/realignment.nf'

workflow {
    REALIGNMENT()
}

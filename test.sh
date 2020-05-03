#!/bin/bash

inputs=./data/*.json

frames=( )

exebase=maph
outdir=./data
expectedoutdir=./data/expected-output
outputext=png
use_stdin="false"
use_pushpop="false"
#use_defaultgen="true"

#===============================================================================

source ./submodules/bat/test.sh


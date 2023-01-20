#!/bin/bash

function run_a1_at_a2QAM_a3xa4_output_a5() {
    $1 $2 1000 $3 $4 1 | grep "SimResult" >> $5-noise.res.csv
    $1 $2 1000 $3 $4 0 | grep "SimResult" >> $5.res.csv
    $1 $2 1000 $3 $4 1 | grep "SimResult" >> $5-noise.res.csv
    $1 $2 1000 $3 $4 0 | grep "SimResult" >> $5.res.csv
    $1 $2 1000 $3 $4 1 | grep "SimResult" >> $5-noise.res.csv
    $1 $2 1000 $3 $4 0 | grep "SimResult" >> $5.res.csv
    $1 $2 1000 $3 $4 1 | grep "SimResult" >> $5-noise.res.csv
    $1 $2 1000 $3 $4 0 | grep "SimResult" >> $5.res.csv
    $1 $2 1000 $3 $4 1 | grep "SimResult" >> $5-noise.res.csv
    $1 $2 1000 $3 $4 0 | grep "SimResult" >> $5.res.csv
}

QAM=("4" "16" "64")
SxR=("4" "6" "8" "12")

function batch_test_a1exe_output_a2() {
    echo "tag,QAM,sender,receiver,sim_time,SymbolER,SendER,exe_time" >> $2.res.csv
    echo "tag,QAM,sender,receiver,sim_time,SymbolER,SendER,exe_time" >> $2-noise.res.csv
    for qam in ${QAM[@]}; do
        for sxr in ${SxR[@]}; do
            run_a1_at_a2QAM_a3xa4_output_a5 $1 $qam $sxr $sxr $2
        done
    done
}

rm -f *.res.csv

# Zero Forcing
make clean
make OPT+=-DZF
batch_test_a1exe_output_a2 ./mysim.exe zf

# Sphere Decoding
make clean
make OPT+=-DSP
batch_test_a1exe_output_a2 ./mysim.exe sp

# Sphere Decoding with Radius Initializing Optimization
make clean
make OPT+=-DSP OPT+=-DSP_RADIUS_OPT
batch_test_a1exe_output_a2 ./mysim.exe sp-init-opt
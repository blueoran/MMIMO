#!/bin/bash

function run_exe_qam_send_recv_output_SNR() {
    $1 $2 1000 $3 $4 1 $6 | grep "SimResult" >> $5-noise.res.csv
    $1 $2 1000 $3 $4 1 $6 | grep "SimResult" >> $5-noise.res.csv
    $1 $2 1000 $3 $4 1 $6 | grep "SimResult" >> $5-noise.res.csv
    $1 $2 1000 $3 $4 1 $6 | grep "SimResult" >> $5-noise.res.csv
    $1 $2 1000 $3 $4 0 $6 | grep "SimResult" >> $5.res.csv
}

QAM=("4" "16" "64")
SxR=("4" "6" "8")

function batchtest_exe_output() {
    echo "tag,QAM,sender,receiver,sim_time,SymbolER,SendER,exe_time" >> $2.res.csv
    echo "tag,QAM,sender,receiver,SNR(dB),S/N,sim_time,SymbolER,SendER,exe_time" >> $2-noise.res.csv
    for qam in ${QAM[@]}; do
        for sxr in ${SxR[@]}; do
            run_exe_qam_send_recv_output_SNR $1 $qam $sxr $sxr $2 -0.3
            run_exe_qam_send_recv_output_SNR $1 $qam $sxr $sxr $2 -0.2
            run_exe_qam_send_recv_output_SNR $1 $qam $sxr $sxr $2 -0.1
            run_exe_qam_send_recv_output_SNR $1 $qam $sxr $sxr $2 0
            run_exe_qam_send_recv_output_SNR $1 $qam $sxr $sxr $2 0.1
            run_exe_qam_send_recv_output_SNR $1 $qam $sxr $sxr $2 0.2
            run_exe_qam_send_recv_output_SNR $1 $qam $sxr $sxr $2 0.3
            run_exe_qam_send_recv_output_SNR $1 $qam $sxr $sxr $2 0.4
            run_exe_qam_send_recv_output_SNR $1 $qam $sxr $sxr $2 0.5
            run_exe_qam_send_recv_output_SNR $1 $qam $sxr $sxr $2 1
            run_exe_qam_send_recv_output_SNR $1 $qam $sxr $sxr $2 2
        done
    done
}

rm -f *.res.csv

# Zero Forcing
make clean
make OPT+=-DZF
batchtest_exe_output ./mysim.exe zf

# Sphere Decoding
make clean
make OPT+=-DSP
batchtest_exe_output ./mysim.exe sp

# Sphere Decoding with Radius Initializing Optimization
make clean
make OPT+=-DSP OPT+=-DSP_RADIUS_OPT
batchtest_exe_output ./mysim.exe sp-init-opt
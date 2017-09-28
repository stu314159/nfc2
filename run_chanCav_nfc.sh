#!/bin/bash

# arguments
#
# 1 - N_divs 

MAT_FILE=ChanCavityTest.mat

Num_ts=45001
ts_rep_freq=1000
Warmup_ts=0
plot_freq=5000
Re=30000
dt=0.00001
Cs=50
Restart_flag=0

aprun -n 1 ./channel_cavity_geom.py $1


# add arguments...
aprun -n 1 ./genInput.py $MAT_FILE $Num_ts $ts_rep_freq $Warmup_ts \
$plot_freq $Re $dt $Cs $Restart_flag  

module swap PrgEnv-gnu PrgEnv-pgi
module swap cudatoolkit/5.5.51-1.0502.9594.1.1 cudatoolkit
module load craype-accel-nvidia35

aprun -B ./NFC

#module swap PrgEnv-pgi PrgEnv-gnu

aprun -n 1 ./processNFC.py

#!/bin/bash

echo "Starting HGC Calibration script script"
echo "Running hard coded runs"
echo "Analysing all events"
echo "Coincident runs"
echo "showall"
SCRIPTPATH="/home/${USER}/work/Jlab/hallc_replay_lt/CALIBRATION/shms_hgcer_calib/new_calib_code"
cd $SCRIPTPATH
root -l   "$SCRIPTPATH/run_cal.C (1, -1, 1, \"showall\")"
exit 0 

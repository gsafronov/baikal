#!/bin/bash

if [ $1 = "-proc" ]; then
	nProc="$2"
	echo "$nProc"
fi

for ((j=0; j<4; j++)); do
	iChan=$((j+1))
	echo "run calibration of channel $iChan"
#	root -q -b rec/root -b rec/timeCalib_stage1_mc.C\(\"/home/local1/work/baikal/mc/cors-n8x10m-2sec-tr112-thr4.0-1.5.wout\",\"OUT/timeCalibOut.root\",$iChan,kTRUE\) >& calib_log/log_channel_$iChan &;
	if ((iChan%nProc==0)); then
		echo "yoyo"
	fi
done

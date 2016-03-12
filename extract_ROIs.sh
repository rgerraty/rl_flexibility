#!/bin/bash

####CODE for extracting timecourses from ROIs, given 4D image and folder containing ROIs
####Note: 4D image must have previously been registered to standard space
####images should be in .nii.gz format, but this can be easily changed

4d_data=$1
rois=$2

if [ -z $2 ];
	then
	echo You must supply a 4D timeseries AND a folder containing regions to extract
	echo Example:
	echo ~/GitHub/rl_flexibility/extract_ROIs.sh res4d.nii.gz ~/Harvard-Oxford_ROIs/
else
	ts_dir=$(dirname $4d_data)/$(basename $rois)
	mkdir $ts_dir
	for r in $rois/*nii.gz
		do
			fslmeants -i $4d_data -o $ts_dir/$(basename $r .nii.gz).txt -m $r --eig
	done
	past 

fi
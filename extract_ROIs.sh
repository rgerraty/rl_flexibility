#!/bin/bash

####CODE for extracting timecourses from ROIs, given 4D image and folder containing ROIs
####Note: 4D image must have previously been registered to standard space
####images should be in .nii.gz format, but this can be easily changed

fourd_data=$1
rois=$2

if [ -z $2 ];
	then
	echo You must supply a 4D timeseries AND a folder containing regions to extract
	echo Example:
	echo ~/GitHub/rl_flexibility/extract_ROIs.sh res4d.nii.gz ~/Harvard-Oxford_ROIs/
else
	if [ ! -z $3 ]
		then
		ts_dir=$3
	else
		ts_dir=$(dirname $fourd_data)/$(basename $rois)
	fi

	if [ -d $ts_dir ]
		then
		echo $ts_dir already exists!
	else
		mkdir $ts_dir
		for r in $rois/*nii.gz
		do
			fslmeants -i $fourd_data -o $ts_dir/$(basename $r .nii.gz).txt -m $r --eig
		done
		paste $ts_dir/*txt
	fi
fi
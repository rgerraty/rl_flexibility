#!/bin/bash

####CODE for running 1st level confound analysis, need to have successfully run initial preprocessing and created confound regressors to run
if [ -z $1 ]; then
echo Script for running 1st level confound regression\. Make sure to enter which confound files to include in analysis'\n'Usage\:'\n'/home/rgerraty/1st_level_conf.sh 36par+spikes.txt

else
#paths to feat directories containing preprocessed data
paths=/data/engine/kfoerde/FoodChoice/restingstate/*/*/*_237vol_120s_5mm_wm.feat/

for i in $paths
	do
	#for each preprocessing directory check for filtered data, get TR and number of volumes, create and run .fsf file from template 
	if [ -e $i$1 ];then 
	conf=$i$1
	out=$i../`basename $1 .txt`
	filtdata=`ls $i\filtered_func_data.nii.gz`
	TR=`fslinfo $filtdata | grep pixdim4 | awk '{ print $2}'`
	#get TR
	vols=`fslinfo $filtdata | grep ^dim4 | awk '{ print $2}'`
	#get number of volumes
	#replace dummy lines in template fsf to make subject-specific temp fsf file
	sed -e 's:XXOUTPUTXX:'$out':g' -e 's:XXTRXX:'$TR':g' -e 's:XXVOLSXX:'$vols':g' -e 's:XX4DDATAXX:'$filtdata':g' -e 's:XXCONFXX:'$conf':g'</home/rgerraty/scripts/fsf_files/conf_reg_design.fsf>tmp.fsf
	feat tmp.fsf #run temp file
	rm -rf tmp.fsf

	else echo -e $1 does not exist in directory $i.'\n'moving on to next run
	fi
done

fi
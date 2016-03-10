#!/bin/bash

####CODE for running 1st level confound analysis, need to have successfully run initial preprocessing and created confound regressors to run
if [ -z $2 ]; then
echo -e Script for running 1st level confound regression\. 
echo -e Make sure to enter input volumes as well as which confound files to include in analysis'\n'Usage\:
echo -e '\n'~/GitHub/rl_flexibility/1st_level_conf.sh filtered_func_data.nii.gz 36par+spikes.txt

else
#paths to feat directories containing preprocessed data


#check for filtered data, get TR and number of volumes, create and run .fsf file from template 
if [ -e $1 ];
	then 
	if [ -e $2 ];
		then
		conf=$(readlink -e $2)
		filtdata=$(readlink -e $1)
		cp $filtdata $(dirname $filtdata)/../$(basename $filtdata)_tmp
		filtdata_tmp=$(readlink -e $(dirname $filtdata)/../$(basename $filtdata)_tmp)
		out=$(dirname $filtdata)/../$(basename $conf .txt)
		
		#get TR
		TR=`fslinfo $filtdata | grep pixdim4 | awk '{ print $2}'`

		#get number of volumes
		vols=`fslinfo $filtdata | grep ^dim4 | awk '{ print $2}'`
		
		#replace dummy lines in template fsf to make subject-specific temp fsf file
		sed -e 's:XXOUTPUTXX:'$out':g' -e 's:XXTRXX:'$TR':g' -e 's:XXVOLSXX:'$vols':g' -e 's:XX4DDATAXX:'$filtdata_tmp':g' -e 's:XXCONFXX:'$conf':g'<~/GitHub/rl_flexibility/conf_reg_design.fsf>tmp.fsf
		feat tmp.fsf #run temp file
		rm -rf tmp.fsf
		rm -rf $filtdata_tmp
	
	else 
		echo -e $2 does not exist
	
	fi

else 
	echo -e $1 does not exist.
fi

fi
#!/bin/bash


featpath=$1 #path to feat directories containing preprocessed data
structpath=$2 #relative path from feat directory to structural folder with wm, csf images

if [ -z $2 ]
	then
	echo Code for extracting confound timecourses from preprocessed data
	echo Need to provide feat directory as well as anatomical directory with csf and white matter masks
	echo Can also provide z-score cut-off for high-motion timepoints \(spikes\)
	echo Usage:
	echo ~/GitHub/rl_flexibility/fsl_extract_confts.sh Learn1_PEprior.feat/ structural/mprage.anat 2.5
else
	if [ -z $3 ];
		then
		sd_thresh=2.5 #outlier threshold in standard deviations from global mean of relative motion
	else
		sd_thresh=$3;
	fi


	struct=`readlink -e $structpath`;
	feat=`readlink -e $featpath`;

	echo generating confound masks...
	fslmaths `ls $struct/*_pve_2.nii.gz` -thr .75 -bin $struct/wm_thr75.nii.gz;
	fslmaths `ls $struct/*_pve_0.nii.gz` -thr .75 -bin $struct/csf_thr75.nii.gz;
	flirt -in $struct/wm_thr75.nii.gz -ref $feat/reg/example_func.nii.gz -applyxfm -init $feat/reg/highres2example_func.mat -out $feat/reg/wm; 
	flirt -in $struct/csf_thr75.nii.gz -ref $feat/reg/example_func.nii.gz -applyxfm -init $feat/reg/highres2example_func.mat -out $feat/reg/csf;
	fslmeants -i $feat/filtered_func_data.nii.gz -o $feat/csf.txt -m $feat/reg/csf; 
	fslmeants -i $feat/filtered_func_data.nii.gz -o $feat/wm.txt -m $feat/reg/wm;
	fslmeants -i $feat/filtered_func_data.nii.gz -o $feat/wb.txt -m $feat/mask
	paste -d "  " $feat/wb.txt $feat/wm.txt $feat/csf.txt $feat/mc/prefiltered_func_data_mcf.par > $feat/9par.txt;
	~/GitHub/rl_flexibility/mp_diffpow_rtg.sh $feat/9par.txt $feat/9par_sq_der_sqder
	/usr/share/fsl/5.0/bin/mp_diffpow.sh $feat/mc/prefiltered_func_data_mcf.par $feat/6par_sq_der_sqder
	cut -f22-40 -d\  $feat/9par_sq_der_sqder.dat>$feat/9par_diff.txt;
	cut -f28-40 -d\  $feat/9par_sq_der_sqder.dat>$feat/6par_diff.txt;
	paste -d "  " $feat/9par.txt $feat/9par_sq_der_sqder.dat > $feat/36_par_conf.txt;
	paste -d "  " $feat/9par.txt $feat/9par_diff.txt> $feat/18_par_conf.txt;
	paste -d "  "  $feat/mc/prefiltered_func_data_mcf.par $feat/6par_diff.txt> $feat/12_par_conf.txt;
	paste -d "  " $feat/mc/prefiltered_func_data_mcf.par $feat/6par_sq_der_sqder.dat > $feat/24_par_conf.txt;


	motion_mean=`awk '{sum+=$1} END {print sum/NR}' $feat/mc/prefiltered_func_data_mcf_rel_mean.rms`
	motion_SD=`awk '{delta = $1 - avg; avg += delta / NR; mean2 += delta * ($1 - avg); } END { print sqrt(mean2 / NR); }' $feat/mc/prefiltered_func_data_mcf_rel_mean.rms`
	spike_thresh=`echo $motion_mean + $sd_thresh*$motion_SD | bc`
	rel_rms=$feat/mc/prefiltered_func_data_mcf_rel.rms
	paste $rel_rms>tmp.txt
	matlab -nosplash -nojvm -r "addpath ~/GitHub/rl_flexibility/;make_spike_regs('`echo $feat`','tmp.txt',`echo $spike_thresh`, `echo $sd_thresh`);quit"
	rm -rf tmp.txt

	if [ -e $feat/spike_regressors_$sd_thresh\_SD.txt ];
		then 
		paste -d "  " $feat/36_par_conf.txt $feat/spike_regressors_2.5_SD.txt>$feat/36par+spikes.txt; 
		paste -d "  " $feat/18_par_conf.txt $feat/spike_regressors_2.5_SD.txt >$feat/18par+spikes.txt;
		paste -d "  " $feat/24_par_conf.txt $feat/spike_regressors_2.5_SD.txt >$feat/24par+spikes.txt;
		paste -d "  " $feat/12_par_conf.txt $feat/spike_regressors_2.5_SD.txt>$feat/12par+spikes.txt;
	else 
		for j in `ls $feat/*par_conf.txt`; 
			do cp $j $feat/`basename $j | cut -f1 -d_`par+spikes.txt
		done
	fi
fi
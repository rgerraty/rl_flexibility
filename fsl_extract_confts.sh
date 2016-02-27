#!/bin/bash


paths=/data/engine/kfoerde/FoodChoice/restingstate/*/*/*_237vol_120s_5mm_wm.feat/ #paths to feat directories containing preprocessed data
structpath=../structural/.anat #relative path from feat directory to structural folder with wm, csf images
sd_thresh=2.5 #outlier threshold in standard deviations from global mean of relative motion


for i in $paths; 
	do struct=`ls -d $i$structpath`;
	echo $i
	fslmaths `ls $struct/*_pve_2.nii.gz` -thr .75 -bin $struct/wm_thr75.nii.gz;
	fslmaths `ls $struct/*_pve_0.nii.gz` -thr .75 -bin $struct/csf_thr75.nii.gz;
	flirt -in $struct/wm_thr75.nii.gz -ref $i\reg/example_func.nii.gz -applyxfm -init $i\reg/highres2example_func.mat -out $i\reg/wm; 
	flirt -in $struct/csf_thr75.nii.gz -ref $i\reg/example_func.nii.gz -applyxfm -init $i\reg/highres2example_func.mat -out $i\reg/csf;
	fslmeants -i $i\filtered_func_data.nii.gz -o $i\csf.txt -m $i\reg/csf; 
	fslmeants -i $i\filtered_func_data.nii.gz -o $i\wm.txt -m $i\reg/wm;
	fslmeants -i $i\filtered_func_data.nii.gz -o $i\wb.txt -m $i\mask
	paste -d "  " $i\wb.txt $i\wm.txt $i\csf.txt $i\mc/prefiltered_func_data_mcf.par > $i\9par.txt;
	/home/rgerraty/scripts/mp_diffpow_rtg.sh $i\9par.txt $i\9par_sq_der_sqder
	/usr/share/fsl/5.0/bin/mp_diffpow.sh $i\mc/prefiltered_func_data_mcf.par $i\6par_sq_der_sqder
	cut -f22-40 -d\  $i\9par_sq_der_sqder.dat>$i\9par_diff.txt;
	cut -f28-40 -d\  $i\9par_sq_der_sqder.dat>$i\6par_diff.txt;
	paste -d "  " $i\9par.txt $i\9par_sq_der_sqder.dat > $i\36_par_conf.txt;
	paste -d "  " $i\9par.txt $i\9par_diff.txt> $i\18_par_conf.txt;
	paste -d "  "  $i\mc/prefiltered_func_data_mcf.par $i\6par_diff.txt> $i\12_par_conf.txt;
	paste -d "  " $i\mc/prefiltered_func_data_mcf.par $i\6par_sq_der_sqder.dat > $i\24_par_conf.txt;
done

motion_mean=`awk '{sum+=$1} END {print sum/NR}' $paths\mc/prefiltered_func_data_mcf_rel_mean.rms`
motion_SD=`awk '{delta = $1 - avg; avg += delta / NR; mean2 += delta * ($1 - avg); } END { print sqrt(mean2 / NR); }' $paths\mc/prefiltered_func_data_mcf_rel_mean.rms`
spike_thresh=`echo $motion_mean + $sd_thresh*$motion_SD | bc`
rel_rms=$paths\mc/prefiltered_func_data_mcf_rel.rms
paste $rel_rms>tmp.txt
matlab -nosplash -nojvm -r "addpath ~/GitHub/rl_flexibility/;make_spike_regs('`echo $paths`','tmp.txt',`echo $spike_thresh`, `echo $sd_thresh`);quit"
rm -rf tmp.txt

for i in $paths;
do if [ -e $i\spike_regressors_2.5_SD.txt ];
	then paste -d "  " $i\36_par_conf.txt $i\spike_regressors_2.5_SD.txt>$i\36par+spikes.txt; 
		paste -d "  " $i\18_par_conf.txt  $i\spike_regressors_2.5_SD.txt >$i\18par+spikes.txt;
		paste -d "  " $i\24_par_conf.txt $i\spike_regressors_2.5_SD.txt >$i\24par+spikes.txt;
		paste -d "  " $i\12_par_conf.txt $i\spike_regressors_2.5_SD.txt>$i\12par+spikes.txt;
	else 
		for j in `ls $i*par_conf.txt`; 
			do cp $j $i`basename $j | cut -f1 -d_`par+spikes.txt
		done
	fi
done
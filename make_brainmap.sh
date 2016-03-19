val_ind=$1
region_list=$2
val_list=$3
p=$4

if [ -z $3 ];
then
	echo code for generating statistical maps from thresholded ROIs
	echo need to specify list of ROIs, list of statistics, and indices of regions passing threshold
	echo ~/GitHub/rl_flexibility/make_brainmap.sh p005_regions.txt HO_name_list.txt p_vals_learning_glm.csv 
	echo Raphael Gerraty 2016
else

k=1;
for i in $(cat $val_ind); 
	do 

	region=$(cat $region_list | awk -v roi=$i 'FNR==roi { print ;exit }'); 
	val=$(cat $val_list | awk -v val=$i 'FNR==val { print;exit }');
	echo $region $val

	if [ "$p" == "p" ];
	then
		val="$(sed 's/[eE]+\{0,1\}/*10^/g' <<<"$val")"
		val=$(echo 1 - $val | bc -l) 
	fi

	if [ $k -eq 1  ];
	then
		fslmaths $region -mul $val $(basename $val_ind .txt)
	else
		fslmaths $region -mul $val roi_tmp
		fslmaths $(basename $val_ind .txt).nii.gz -add roi_tmp.nii.gz $(basename $val_ind .txt) 
	fi
	k=$(($k+1))
done
rm -rf roi_tmp.nii.gz
fi
	


#Network Flexibility and Reinforcement Learning
##Raphael Gerraty, 2015-2016

Descriptions and example scripts for running network preprocessing and analysis functions contained in this repository. See paper for details when it comes out.



### Extended Preprocessing
Because of the known effect of motion on measures of connectivity, we followed up standard preprocessing in FSL with an extended nuisance regression. Affine transformation parameters from motion correction, CSF, white matter, and whole-brain signals are regressed against preprocessed 4D data, along with the squares, derivatives, and squared derivatives of these confounds. See Satterthwaite et al 2013 for details. 

```{.bash}
for i in /data/engine/rgerraty/learn_dyncon/4*/Learn*/filtered_func_data.nii.gz
  do
  subdir=$(dirname $i)

  #extract confound timecourses from preprocessed data
  #need to provide feat directory as well as anatomical directory
  #can also provide z-score cut-off for high-motion timepoints (spikes)
  ~/GitHub/rl_flexibility/fsl_extract_confts.sh $subdir\
   $subdir/../structural/mprage.anat 3

  #run 1st level confound regression, using template .fsf and confound files
  ~/GitHub/rl_flexibility/1st_level_conf.sh $i $subdir/36par+spikes.txt
done
```




### Running nonlinear registration with FNIRT
After nuisance regression has been run, the residual timeseries needs to be transformed into standard space (in this case, MNI). Make sure fsl\_anat has been run on each structural image first. The following bash code was used to perform these transformations:

``` {.bash}
#fnirt has already been run, just applying transformation
for i in /data/engine/rgerraty/learn_dyncon/4*/Learn*; 
    do 
    #run linear registration on example functional image
    flirt -ref $i/../structural/mprage.anat/T1_biascorr_brain.nii.gz\
     -in $i/reg/example_func.nii.gz\
     -omat $i/reg/example_func2highres.mat;

     echo warping $i;

    #apply warp from FNIRT to preprocessed 4D data
    applywarp --ref=$FSLDIR/data/standard/MNI152_T1_2mm.nii.gz\
      --in=$i/36par+spikes.feat/stats/res4d.nii.gz\
      --out=$i/36par+spikes.feat/stats/res4d_std.nii.gz\
      --warp=$i/../structural/mprage.anat/T1_to_MNI_nonlin_field.nii.gz\
      --premat=$i/reg/example_func2highres.mat;  
done
```



### Extracting time courses
Once the preprocessed images have been registered, we extract mean timecourses for each Harvard-Oxford ROI, using the function extract_ROIs.sh. The output of this function is a timecourse for each ROI in the specified input folder, as well as a .txt file containing all of the ROIs. The bash code used to run this function on each learning block for each subject is below:

``` {.bash}
for i in /data/engine/rgerraty/learn_dyncon/4*/Learn?_PEprior.feat/36par+spikes.feat/; 
    do 
    #extract timeseries (mean or 1st eigenvector, see function) data from each ROI in ~/Harvard-Oxford_ROIs/ 
    ~/GitHub/rl_flexibility/extract_ROIs.sh $i/stats/res4d_std.nii.gz ~/Harvard-Oxford_ROIs/ $i/H-O_rois/;
done
```



### Calculate coherence matrices for each time window
Connectivity between pairs of ROIs was measured by average magnitude squared coherence in the .06-.12 Hz band, computed in MATLAB. The code below calls a function for creating a coherence matrices in a specified frequency range for specified time windows (in this case 25 TRs, or 50 s). These are saved as a .mat file for multi-slice community detection. 


``` {.matlab}
addpath ~/GitHub/rl_flexibility
%read in all subject/run ROI timeseries directories 
[a,b]=system('ls -d /data/engine/rgerraty/learn_dyncon/4*/Learn?_PEprior.feat/36par+spikes.feat/H-O_rois');
c=strread(b,'%s');

for i=1:size(c,1)

    %calculate coherence per time window from concatenated ROI file
    filename=char(strcat(c(i),'/all_rois.txt'))
    %need to specify filename, window length in TR, sampling rate, bandpass 
    conn_cell=coherence_by_block(filename,25,.5,.06,.12);
    save(char(strcat(c(i),'/conn_cells')),'conn_cell')

end
```




### Run multi-slice community detection and flexibility statistics

Input coherence matrix for each block. Also need number of blocks,
resolution and coupling parameters. In Matlab

``` {.matlab}
%need multi-slice, flexibility codes not yet on GitHub for network_diags to run 
addpath ~/GitHub/rl_flexibility
addpath ~/scripts/MATLAB/GenLouvain_for_Raphael/
addpath ~/scripts/MATLAB/Bassett_Code/

%read in data
[a,b]=system('ls -d /data/engine/rgerraty/learn_dyncon/4*/Learn?_PEprior.feat/36par+spikes.feat/H-O_rois/');
c=strread(b,'%s');

%concatenate runs for each subject
numruns=4
k=1;
for j=1:size(c,1)/numruns
    c(k)
    conn_cell_cat=[];
    for i=1:numruns 
        load(strcat(char(c(k-1+i)),'/conn_cells'))
        conn_cell_cat=cat(3,conn_cell_cat,conn_cell)
    end

    %network_diags code:
    %runs multi-slice community detection
    %gives flexibility for each run
    %also allegiance matrix (not using yet)
    %need to specify number of blocks, simulations, coupling, resolution
    [a_mat,flex]=network_diags(conn_cell_cat,4,500,1,1.1813)
    save(char(strcat(c(k),'/../../../a_mat')),'a_mat')
    save(char(strcat(c(k),'/../../../flex')),'flex')
    k=k+numruns;
end
```



### Pull flexibility statistics

For plotting and preparing for heirarchical models. Matlab.

``` {.matlab}

%load data and concatenate flexibility statistics
[a,b]=system('ls -d /data/engine/rgerraty/learn_dyncon/4*/flex.mat');
c=strread(b,'%s');
flex_cat=[];
for j=1:size(c,1)
    load(char(c(j)))
    flex_cat=cat(3,flex_cat,flex)
end
plot(squeeze(mean(flex_cat)))

block=repmat([1:4]',22,1);
sub=repmat([1:22]',1,4)'
sub=sub(:);

%reshape whole-brain average flexibility
meanflex=squeeze(mean(flex_cat));
meanflex=meanflex(:);

%get striatal average flexibility
%check to make sure indices are correct
[trash,roi_names]=system('ls  ~/Harvard-Oxford_ROIs/*nii.gz | xargs -n1 basename');
roi_names=strread(roi_names,'%s');
str_ind=[49,51,54,104,106,109];
roi_names(str_ind)

strflex=squeeze(mean(flex_cat(str_ind,:,:)));
strflex=strflex(:);

plot(squeeze(mean(flex_cat(str_ind,:,:))))

%write out csv for modeling in R
flexdata=[sub block meanflex strflex]
dlmwrite('/data/engine/rgerraty/learn_dyncon/flexdata.csv',flexdata) 

%get flexibility scores for each ROI for each run for whole-brain search
%can prob do this more effeciently but this is easier to see, harder to botch
flex_allrois=[];
for i=1:size(flex_cat,3)
  flex_allrois=[flex_allrois;flex_cat(:,:,i)'];
end

dlmwrite('/data/engine/rgerraty/learn_dyncon/flex_allrois.csv',flex_allrois) 
```

REML and fully Bayesian models for the effect of striatal and whole-brain flexibility on reinforcement learning 

```{.r}
library(reshape2)
library(lme4)
library(rstanarm)
library(brms)

#prepare data for binomial logistic modelling

#read in trial-by-by trial behavioral data 
data<-read.delim('/data/engine/rgerraty/learn_dyncon/behav_data.tsv',header=1)
data$block<-rep(rep(seq(1,4,1),each=30),22)

#read in block-level flexibility data
flexdat<-read.csv('/data/engine/rgerraty/learn_dyncon/flexdata500sim.csv',header=0)
names(flexdat)<-c('subject','block','wb_flex','str_flex')

#get proportion correct for each block
flex_behav$correct<- melt(tapply(data$optCor,list(data$block,data$subjectNum),mean,na.rm=1))$value

#get weights (number of trials with responses) for each block
nacount<-is.na(data$optCor)
weights<-melt(tapply(nacount,list(data$block,data$subjectNum),sum))
flex_behav$weights<-30-weights[,3]

#write out for future use
flex_behav<-cbind(flex_behav,flexdat)
write.csv(flex_behav,'/data/engine/rgerraty/learn_dyncon/flex_behav.csv')


#fit mixed-effects model in lme4 with wald approximation to p-value
mlearn_glmer<-glmer(correct~str_flex+(str_flex||subject),data=flex_behav,family=binomial,weights=flex_behav$weights)

#for posterior inference, run bayesian model using brms wrapper for stan
flex_behav$numcorr<-flex_behav$correct*flex_behav$weights
mlearn_stan<-brm(correct~str_flex+(str_flex|subject),data=flex_behav,family=binomial)


```


Whole-brain search for effects of flexibility on learning


```{.r}
#load in data
library(lme4)
roi_data<-read.csv('/data/engine/rgerraty/learn_dyncon/flex_allrois.csv',header=0)
behav<-read.csv('/data/engine/rgerraty/learn_dyncon/flex_behav.csv')

p<-NULL
b<-NULL


for(i in 1:dim(roi_data)[2]){
  mtmp<-summary(glmer(behav$correct~roi_data[,i]+(roi_data[,i]||behav$subject),family=binomial,weights=behav$weights))
  if(dim(summary(mtmp)[10]$coefficients)[2]>=4){
    p[i]<-summary(mtmp)[10]$coefficients[2,4]
    b[i]<-summary(mtmp)$coefficients[2,1]
  } else{
    p[i]<-NA
    b[i]<-NA
  }
}
write(p,'/data/engine/rgerraty/learn_dyncon/p_vals_learning_glm.csv')
write(b,'/data/engine/rgerraty/learn_dyncon/b_vals_learning_glm.csv')
```

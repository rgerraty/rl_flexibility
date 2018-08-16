# Reinforcement learning model

### Combine raw reinforcement learning data from individual subjects
Matlab
```.matlab
%load in the data
a=dir('/data/engine/rgerraty/learn_dyncon/behavior/*_tb_l*_svlo*.mat');

a={a.name};

%choice - 1 & 0 for right and left flower "resp"
%payoff - 0 incorr, 1 corr "shown_corr"

longform=[];
for k = 1:length(a)
	k
	clear subnum chose_right shown_corr trial
	filea=a{(k)};
	load(a{(k)});	
    subnum=str2double(num2str(a{k}(1:3)))';
    chose_right=double(strcmp(resp,'b'))';
    shown_corr=double(shown_corr)';
    no_resp=shown_corr==2;
    chose_right(no_resp)=NaN;
    shown_corr(no_resp)=NaN;
    ntrials=size(resp',1);
    trial=1:ntrials';
    longform=[longform;repmat(subnum,ntrials,1) stim_shown' trial' chose_right shown_corr];
end
header={'sub','stim','trial','choice','fb'};
header2=sprintf('%s,',header{:});header2(end)=[];
dlmwrite('/data/engine/rgerraty/learn_dyncon/behavior/choice_fb_long.csv',...
	header2,'')

dlmwrite('/data/engine/rgerraty/learn_dyncon/behavior/choice_fb_long.csv',...
	longform,'-append','delimiter',',')

!cp /data/engine/rgerraty/learn_dyncon/behavior/choice_fb_long.csv ~/GitHub/rl_flexibility/RL_model/choice_fb_long.csv
```

### Prepare data for stan and run hierarchical bayesian model
R
```.r
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

dat <- read.csv('~/GitHub/rl_flexibility/RL_model/choice_fb_long.csv')

choices<-unique(na.omit(dat$choice))

dat$chosen<-dat$choice
dat$chosen[dat$choice==choices[1]]=1
dat$chosen[dat$choice==choices[2]]=2
dat$unchosen=abs(dat$chosen-3)
dat$chose_two<-dat$chosen==2

NR=4;
subs = unique(dat$sub);
NS = length(subs);
NStim=length(unique(dat$stim));
MT=max(dat$trial);
NT = array(0,NS);
choice = array(0,c(NS,MT));
unchoice=choice;
choice_two=choice;
rew = choice;
stim=choice;
run=choice;

dat$run<-ceiling(dat$trial/MT*NR)

for (i in 1:NS) {
  NT[i] = nrow(subset(dat,sub==subs[i]));

  stim[i,1:NT[i]] = subset(dat,sub==subs[i])$stim;
  run[i,1:NT[i]] = subset(dat,sub==subs[i])$run;
  
  #choice and reward history
  choice[i,1:NT[i]] = subset(dat,sub==subs[i])$chosen;
  unchoice[i,1:NT[i]] = subset(dat,sub==subs[i])$unchosen;

  rew[i,1:NT[i]] = subset(dat,sub==subs[i])$fb;

  #based on choosing second option
  choice_two[i,1:NT[i]] = subset(dat,sub==subs[i])$chose_two;
}

choice[is.na(choice)]<--1
unchoice[is.na(unchoice)]<--1
rew[is.na(rew)]<--1
choice_two[is.na(choice_two)]<--1
rew[run==0]<- -1 

rl_standata = list(NS=NS, NC=2, NStim=NStim, MT=MT, NT= NT, choice=choice, 
stim=stim,choice_two=choice_two,rew=rew );
rl_fit <- stan(file = '~/GitHub/rl_flexibility/RL_model/standard_rl.stan', 
data = rl_standata, iter = 1250, warmup = 250, chains = 4)
save(rl_fit,file='~/Documents/RL_Flexibility/data/rl_fit')

rl_standata_multialpha = list(NS=NS, NC=2, NStim=NStim, MT=MT, NT= NT, choice=choice, 
stim=stim,choice_two=choice_two,rew=rew,run=run,NR=NR );
rl_multialpha_fit <- stan(file = '~/GitHub/rl_flexibility/RL_model/multialpha_rl.stan', 
data = rl_standata_multialpha, iter = 1250, warmup = 250, chains = 4)
save(rl_multialpha_fit,file='~/Documents/RL_Flexibility/data/rl_multialpha_fit')

library(loo)
rl_fit_llik<-extract_log_lik(rl_fit)
rl_multialpha_fit_llik<-extract_log_lik(rl_multialpha_fit)

rl_loo<-loo(rl_fit_llik)
rl_multialpha_loo<-loo(rl_multialpha_fit_llik)

compare(rl_loo,rl_multialpha_loo)
```

### Extract parameters from model
R
```.r
fit_rl<-load('~/Documents/rl_fit')
fit_extract<-extract(rl_fit,permute=T)

betas<-apply(fit_extract$beta,2,mean)
alphas<-apply(fit_extract$alpha,2,mean)

fit_extract$beta_mean<-fit_extract$b1*fit_extract$b2
fit_extract$alpha_mean<-fit_extract$a1/(fit_extract$a1+fit_extract$a2)

str_acorrel<-matrix(0,4000,6)
str_bcorrel<-matrix(0,4000,6)
str_apbcorrel<-matrix(0,4000,6)
str_bpacorrel<-matrix(0,4000,6)
wb_acorrel<-NULL
wb_bcorrel<-NULL
wb_apbcorrel<-NULL
wb_bpacorrl<-NULL


flexdat<-read.csv("~/Documents/flex_allrois.csv",header=F)
flex_behav<-read.csv("~/Documents/flex_behav.csv",header=T)
flexdat$Sub<-rep(seq(1,22,1),each=4)
flexdat$Block<-rep(seq(1,4,1),times=22)
flexdat$Corr<-flex_behav$correct
flexdat$weights<-flex_behav$weights
flexdat<-melt(flexdat,id=c("Sub","Block","Corr","weights"))
names(flexdat)[c(5,6)]<-c("ROI","flex")
flexdat$ROI<-as.factor(as.numeric(flexdat$ROI))
meanflex_rois<-tapply(flexdat$flex,list(flexdat$Sub,flexdat$ROI),mean)
meanflex<-rowMeans(meanflex_rois)

str_ind=c(49,51,54,104,106,109);

for( i in seq(1,dim(fit_extract$beta)[1],1)){
  
  betas_tmp=fit_extract$beta[i,]
  alphas_tmp=fit_extract$alpha[i,]
  k<-1
  
  wb_acorrel[i]<-cor(meanflex,alphas_tmp)
  wb_bcorrel[i]<-cor(meanflex,betas_tmp)
  wb_apbcorrel[i]<-pcor.test(meanflex,alphas_tmp,betas_tmp)[1]
  wb_bpacorrel[i]<-pcor.test(meanflex,betas_tmp,alphas_tmp)[1]
  
  for(j in str_ind){
  str_acorrel[i,k]<-cor(meanflex_rois[,k],alphas_tmp)
  str_bcorrel[i,k]<-cor(meanflex_rois[,k],betas_tmp)
  
  str_bpacorrel[i,k]<-pcor.test(meanflex_rois[,k],betas_tmp,alphas_tmp)$estimate
  str_apbcorrel[i,k]<-pcor.test(meanflex_rois[,k],alphas_tmp,betas_tmp)$estimate
  k<-k+1
  }
}

par(mfrow=c(1,2))

plot(rowMeans(str_bcorrel),rowMeans(str_acorrel),
     col=rgb(0,0,0,alpha=0),pch=21,bg=rgb(0,0,0,alpha=.03),
     xlab="Correlation with Beta",ylab="Correlation with Alpha",
     main="Striatum Flexibility",xlim=c(-.1,.5))
abline(v=0,lty=2)
abline(h=0,lty=2)

plot(wb_bcorrel,wb_acorrel,
     col=rgb(0,0,0,alpha=0),pch=21,bg=rgb(0,0,0,alpha=.03),
     xlab="Correlation with Beta",ylab="Correlation with Alpha",
     main="Whole-brain Flexibility",xlim=c(-.1,.5))
abline(v=0,lty=2)
abline(h=0,lty=2)


Qvals<-apply(fit_extract$Q,c(2,3,4,5),mean)

Q_chosen<-matrix(0,dim(Qvals)[2],dim(Qvals)[1])
Q_unchosen<-Q_chosen

for(i in 1:dim(Qvals)[1]){
  for(j in 1:dim(Qvals)[2]){
      Q_chosen[j,i]<-Qvals[i,j,stim[i,j],choice[i,j]]
      Q_unchosen[j,i]<-Qvals[i,j,stim[i,j],unchoice[i,j]]
  }

}
pe_hyb<-t(apply(fit_extract$delta,c(2,3),mean))



```
###Combine raw reinforcement learning data from individual subjects
```{.matlab}
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
    longform=[longform;repmat(subnum,ntrials,1) trial' chose_right shown_corr];
end
header={'sub','trial','choice','fb'};
header2=sprintf('%s,',header{:});header2(end)=[];
dlmwrite('/data/engine/rgerraty/learn_dyncon/behavior/choice_fb_long.csv',...
	header2,'')

dlmwrite('/data/engine/rgerraty/learn_dyncon/behavior/choice_fb_long.csv',...
	longform,'-append','delimiter',',')

!cp /data/engine/rgerraty/learn_dyncon/behavior/choice_fb_long.csv ~/GitHub/rl_flexibility/RL_model/choice_fb_long.csv
```

###Prepare data for stan and run hierarchical bayesian model
```{.R}
library(rstan)
dat <- read.csv('~/GitHub/rl_flexibility/RL_model/choice_fb_long.csv')

choices<-unique(na.omit(dat$choice))

dat$chosen<-dat$choice
dat$chosen[dat$choice==choices[1]]=1
dat$chosen[dat$choice==choices[2]]=2
dat$unchosen=abs(dat$chosen-3)
dat$chose_two<-dat$chosen==2

subs = unique(dat$sub);
NS = length(subs);
MT=max(dat$trial);
NT = array(0,NS);
choice = array(0,c(NS,MT));
unchoice=choice;
choice_two=choice;
rew = array(0.0,c(NS,MT));

for (i in 1:NS) {
  NT[i] = nrow(subset(dat,sub==subs[i]));
  
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

rl_standata = list(NS=NS, NC=2, MT=MT, NT= NT, choice=choice, 
choice_two=red_choice,rew=rew )
rl_fit <- stan(file = '~/GitHub/hybrid_reinforcement_learning/standard_rl.stan', data = rl_standata, iter = 1250, warmup = 250, chains = 4)
save(standard_fit,file='~/Documents/Hybrid_RL/stanfit_rl')
```
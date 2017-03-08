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
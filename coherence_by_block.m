block_size=25;%size of each block for network analyses
Fs=.5;

base='/data/engine/rgerraty/learn_dyncon/'
rois='/all_HO_rois.txt';

subs=dir(strcat(base,'4*'));

minhz=.06
maxhz=.12

for sub=1:size(subs,1)
	sub
	learn_runs=dir(strcat(base,subs(sub).name,'/','Learn*'));
	num_runs=size(learn_runs,1);
	for run=1:num_runs
		ts_length=size(tss,1);
		k=1;
		for i=1:size(tss,1)/block_size
			conn_mat(:,:,i+(run-1)*size(tss,1)/block_size)=...
			mul_coher(tss(k:k+block_size-1,:),Fs,minhz,maxhz);
			k=k+block_size;
		end
	end
	conn_cells{sub}=mat2cell(conn_mat,size(tss,2),size(tss,2),[ones(1,size(conn_mat,3))]);
end

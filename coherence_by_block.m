function conn_cells=coherence_by_block(filename,block_size,Fs,minhz,maxhz)

%numbers for flexibility paper
%block_size=25;%size of each block for network analyses
%Fs=.5;
%minhz=.06
%maxhz=.12


tss=dlmread(filename);  
ts_length=size(tss,1);
k=1;
for i=1:size(tss,1)/block_size
	conn_mat(:,:,i)=...
	mul_coher(tss(k:k+block_size-1,:),Fs,minhz,maxhz);
	k=k+block_size;
end

conn_cells=mat2cell(conn_mat,size(tss,2),size(tss,2),[ones(1,size(conn_mat,3))]);
end

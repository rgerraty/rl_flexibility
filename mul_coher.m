function coh_mat=mul_coher(tss,Fs,minhz,maxhz)

h=1;
coh_mat=zeros(size(tss,2));

for f=1:size(tss,2)
	for g=h:size(tss,2)
		[coh,freq]=mscohere(tss(:,f),tss(:,g),[],[],128,Fs);
		coh_mat(f,g)=mean(coh(freq<=maxhz & freq>=minhz));
	end
	h=h+1;
end

coh_mat=coh_mat+coh_mat'-eye(size(coh_mat));
end
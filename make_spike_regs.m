
function make_spike_regs(paths, rel_rms_file, rel_thresh, sd_thresh)

paths = strread(paths,'%s','delimiter',' ');
rel_rms=load(rel_rms_file);
for i=1:size(rel_rms,2)
	ind=find(rel_rms(:,i)>rel_thresh); 
	for j=1:size(ind,1)                                                
		spikes(:,j)=zeros(size(rel_rms,1),1);                              
		spikes(ind(j),j)=1;                     
	end
	if size(ind,1)>0
		spikes=[zeros(1,size(spikes,2));spikes];
		dlmwrite(char(strcat(paths(i),'spike_regressors_',num2str(sd_thresh),'_SD.txt')), spikes, 'delimiter', ' ');
		clear spikes
	end
end



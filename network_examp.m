%example for running network diags and workign with output

for s=1:22
	for r=1:4
		[a_mat,flex,S_tmp,Q_tmp]=network_diags(corr_cells{s}(blocks==r),100,mean(best_omega),mean(best_res));
		save(strcat('sub',num2str(s),'_run',num2str(r)),'a_mat','flex')
	end
end

for s=1:22
	for r=1:4
		load(strcat('sub',num2str(s),'_run',num2str(r)))
		flexdat(:,r,s)=flex;
		[i_s(:,r,s),norm_is(:,r,s)]=inter_strength(a_mat,net_lab,3,1);
	end
end

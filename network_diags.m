function [a_mat,flex,S_tmp,Q_tmp]=network_diags(conn_cells,blocks,sim,omega,res)

%Flexibility paper numbers 
%blocks=4
%sim=100
%res=1.1813
%omega=1

for i=1:sim
	[S_tmp(:,:,i), Q_tmp(i)]=multiord_res_norm(conn_cells,omega, res);
	for b=1:blocks
		flex_tmp(:,b,i)=flexibility(S_tmp(:,:,i)');
	end
end;
flex=mean(flex_tmp,3);

for h=1:size(conn_cells,3)
	for i=1:sim
		a_mat_tmp(:,:,i)=zeros(size(S_tmp(:,:,i),1));
		for j=1:size(S_tmp(:,h,i),1)
			for k=1:j
				a_mat_tmp(j,k,i)=S_tmp(j,h,i)==S_tmp(k,h,i);
			end
		end
		a_mat_tmp(:,:,i)=a_mat_tmp(:,:,i)+tril(a_mat_tmp(:,:,i),-1)';
	end
	a_mat(:,:,h)=sum(a_mat_tmp,3)/sim;
end



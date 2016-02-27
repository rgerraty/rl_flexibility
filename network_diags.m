function [a_mat,flex,S_tmp,Q_tmp]=network_diags(conn_cells,sim,omega,res)

for i=1:sim
	[S_tmp(:,:,i), Q_tmp(i)]=multiord_res_norm(conn_cells,omega, res);
	flex_tmp(:,i)=flexibility(S_tmp(:,:,i)');
end;
flex=mean(flex_tmp,2);

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



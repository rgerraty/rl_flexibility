function [q_tmp,q_null]=Q_grid(conn_cell,nsim,res_pars,gamma_pars)

q_tmp=zeros(size(res_pars,1),size(gamma_pars,1),nsim);
q_null=zeros(size(res_pars,1),size(gamma_pars,1),nsim);
	for s=1:nsim
	%grid search paramater optimization based on Q-Qnull difference
			r=1;
			for res = res_pars
			o=1;
  			for g = gamma_pars
    			[~,q_tmp(r,o,s)]=multiord_res_norm(conn_cell,g, res);
    			[~,q_null(r,o,s)]=multiord_res_norm_temporal(conn_cell,g, res);
    			o=o+1;
  			end
  			r=r+1;
  		end
  	end
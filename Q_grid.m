function [q,q_tmp]=Q_grid(conn_cell,nsim,res_pars,c_pars)

q=zeros(size(res_pars,1),size(gamma_pars,1),nsim);
q_tmp=zeros(size(res_pars,1),size(gamma_pars,1),nsim);
	for s=1:nsim
	%grid search paramater optimization based on Q-Qnull difference
			r=1;
			for res = res_pars
			o=1;
  			for c = c_pars
    			[~,q(r,o,s)]=multiord_res_norm(conn_cell,c, res);
    			[~,q_tmp(r,o,s)]=multiord_res_norm_temporal(conn_cell,c, res);
    			o=o+1;
  			end
  			r=r+1;
  		end
  	end
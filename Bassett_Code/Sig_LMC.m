function persistence = Sig_LMC(C, A)
%%
% This a function that using Lumped Markov chain to calculate 
% the significance of clusters in a give communinity structure.
% refer to "Piccardi 2011 in PloS one". 
% Here we normalize the original definition of persistence by 
% the size of the corresponding cluster to get a better 
% INPUT:
% "A" is a N-by-N weighted adjacency matrix
% "C" is a N-by-1 partition(cluster) vector
% OUTPUT:
% normalized persistence probability of all clusters
%%
P = diag(sum(A,2))\A;               % transition matrix
[evec, eval] = eigs(P',1);          
if min(evec)<0
    evec = -evec;
end
pi = evec';                         % stationary distribution
num_node = size(A,1);          
cl_label = unique(C);               % # of clusters
num_cl = numel(cl_label);
H = zeros(num_node, num_cl);
for i = 1:num_cl
    H(:, i) = C==cl_label(i);       % cluster label matrix
end
Q = (diag(pi*H))\H'*diag(pi)*P*H;   % transition matrix of the lumped Markov chain
persistence = diag(Q)./sum(H)'*sum(H(:));

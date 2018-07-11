function P = promiscuity(S)
%PROMISCUITY      Promiscuity coefficient
%
%   P = PROMISCUITY(S) calculates the promiscuity coefficient. The
%   promiscuity of a temporal or multislice network corresponds to the
%   fraction of all the communities in the network in which a node
%   participates at least once.
%
%   Inputs:     S,      pxn matrix of community assignments where p is the
%                       number of slices/layers and n the number of nodes
%
%   Outputs:    P,      Promiscuity coefficient
%
%   Other m-files required: none
%   Subfunctions: none
%   MAT-files required: none
%   _______________________________________________
%   Marcelo G Mattar (08/21/2014) 

[~, numNodes] = size(S);
numCommunities = length(unique(S));
P = zeros(numNodes,1);

for n=1:numNodes
    P(n,1) = (length(unique(S(:,n)))-1) / (numCommunities-1);
	% Notice that P=0 if it only participates in one community and P=1 if
	% it participates in every community of the network
end
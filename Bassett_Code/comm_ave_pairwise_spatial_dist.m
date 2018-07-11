function [ comm_ave_pairwise_spatial_dist_array ] = comm_ave_pairwise_spatial_dist(partitions,locations)
% function to calcuate the average pairwise spatial distance between nodes
% within detected communities

% Written by Sarah Feldt Muldoon
% 
% Output: 
%         comm_ave_pairwise_spatial_dist_array   --- an (M+1)x2 array where M is the number
%         of communities.  The first column contains the community index and
%         the second column contains the average pairwise spatial distance 
%         within that community.  The M+1th entry contains the average pairwise 
%         spatial distance of the network (denoted as community 0).
%  
%
% Input: 
%         'partitions'    --- an Nx1 array where N is the number of nodes
%         in the network.  Each entry contains the communitiy index for
%         node i (where communities are sequentially indexed as 1:M and M
%         is the total number of detected communities).
% 
%         'locations'   ---  a Nxdim array where N is the number of nodes and
%         dim is the spatial dimensionality of the network (1,2,or 3).  The
%         columns contain the (x,y,z) coordinates of the nodes in euclidean
%         space


number_nodes=length(partitions);
number_communities=max(partitions);
comm_ave_pairwise_spatial_dist_array=zeros(number_communities+1,2);
comm_ave_pairwise_spatial_dist_array(:,1)=[1:number_communities 0];

%calcuate for each community
for i=1:number_communities
    community_i_nodes=find(partitions==i);
    number_nodes_in_i=length(community_i_nodes);
    dist_array=[];
    if number_nodes_in_i >= 2
        for j=1:number_nodes_in_i-1
            node_j=community_i_nodes(j);
            for k=j+1:number_nodes_in_i
                node_k=community_i_nodes(k);
                coordinates_j=locations(node_j,:);
                coordinates_k=locations(node_k,:);
                dist_jk=sqrt(sum((coordinates_j-coordinates_k).^2));
                dist_array=[dist_array dist_jk];
            end
        end
        ave_dist_i=mean(dist_array);
        comm_ave_pairwise_spatial_dist_array(i,2)=ave_dist_i;
    end
end

%calculate for the network
dist_array=[];
for i=1:number_nodes-1
    coordinates_i=locations(i,:);
    for j=i+1:number_nodes
        coordinates_j=locations(j,:);
        dist_ij=sqrt(sum((coordinates_i-coordinates_j).^2));
        dist_array=[dist_array dist_ij];
    end
end
ave_dist_network=mean(dist_array);
comm_ave_pairwise_spatial_dist_array(number_communities+1,2)=ave_dist_network;
end



function [ comm_radius_array ] = comm_radius(partitions,locations)
% function to calcuate the radius 
% of detected communities

% Written by Sarah Feldt Muldoon
% 
% Output: 
%         comm_radius_array   --- an (M+1)x2 array where M is the number
%         of communities.  The first column contains the community index and
%         the second column contains radius of that community. 
%         The M+1th entry contains the average community radius of the
%         network (denoted as community 0).
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
comm_radius_array=zeros(number_communities+1,2);
comm_radius_array(:,1)=[1:number_communities 0];
number_nodes_in_communities=zeros(number_communities,1);

%calcuate the radius for each community
for i=1:number_communities
    community_i_nodes=find(partitions==i);
    number_nodes_in_i=length(community_i_nodes);
    number_nodes_in_communities(i)=number_nodes_in_i;
    if number_nodes_in_i >= 2
        position_vectors_i=locations(community_i_nodes,:);
        std_position_vectors_i=std(position_vectors_i,1);
        radius_i=sqrt(sum(std_position_vectors_i.^2,2));
        comm_radius_array(i,2)=radius_i;
    end
end


%calcuate the averge community radius for the partition
std_all_position_vectors=std(locations,1);
radius_all=sqrt(sum(std_all_position_vectors.^2));
average_community_radius=1/(number_nodes*radius_all)*sum(number_nodes_in_communities.*comm_radius_array(1:number_communities,2));
comm_radius_array(number_communities+1,2)=average_community_radius;


end


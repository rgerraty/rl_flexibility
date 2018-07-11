function [ comm_spatial_extent_array ] = comm_spatial_extent(partitions,locations)
% function to calcuate the spatial extent 
% of detected communities

% Written by Sarah Feldt Muldoon
% 
% Output: 
%         comm_spatial_extent_array   --- an (M+1)x2 array where M is the number
%         of communities.  The first column contains the community index and
%         the second column contains the spatial extent 
%         of that community.  The M+1th entry contains the 
%         spatial extent of the network (denoted as community 0). 
%         
%  
%
% Input: 
%         'partitions'    --- an Nx1 array where N is the number of nodes
%         in the network.  Each entry contains the communitiy index for
%         node i (where communities are sequentially indexed as 1:M and M
%         is the total number of detected communities).
% 
%         'locations'   ---  a Nxdim array where N is the number of nodes and
%         dim is the spatial dimensionality of the network (2,or 3).  The
%         columns contain the (x,y,z) coordinates of the nodes in euclidean
%         space
% 



number_nodes=length(partitions);
number_communities=max(partitions);
comm_spatial_extent_array=zeros(number_communities,2);
comm_spatial_extent_array(:,1)=1:number_communities;
number_nodes_in_communities=zeros(number_communities,1);

size_locations=size(locations);
dim=size_locations(2);


%calcuate the spatial extent for each community
for i=1:number_communities
    community_i_nodes=find(partitions==i);
    number_nodes_in_i=length(community_i_nodes);
    number_nodes_in_communities(i)=number_nodes_in_i;
    if dim==2
        [conv_hull_i,area_i]=convhull(locations(community_i_nodes,1),locations(community_i,2));
        spatial_extent_i=area_i/number_nodes_in_i;
    end
    if dim==3
        [conv_hull_i,volume_i]=convhull(locations(community_i_nodes,1),locations(community_i_nodes,2),locations(community_i_nodes,3));
        spatial_extent_i=volume_i/number_nodes_in_i;
    end
    comm_spatial_extent_array(i,2)=spatial_extent_i;
    
end

%calcuate the spatial extent for the entire network

if dim==2
    [conv_hull_i,area_i]=convhull(locations(:,1),locations(:,2));
    spatial_extent_network=area_i/number_nodes;
end
if dim==3
    [conv_hull_i,volume_i]=convhull(locations(:,1),locations(:,2),locations(:,3));
    spatial_extent_network=volume_i/number_nodes;
end
comm_spatial_extent_array(number_communities+1,2)=spatial_extent_network;


end


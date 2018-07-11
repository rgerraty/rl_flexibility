function [ comm_laterality_array ] = comm_laterality(partitions,categories)
% function to calcuate the laterality 
% of detected communities

% Written by Sarah Feldt Muldoon
% 
% Output: 
%         comm_laterality_array   --- an (M+1)x2 array where M is the number
%         of communities.  The first column contains the community index and
%         the second column contains the laterality of that community. 
%         The M+1th entry contains the  laterality of the partition of the
%         network (denoted as community 0).
%  
%
% Input: 
%         'partitions'    --- an Nx1 array where N is the number of nodes
%         in the network.  Each entry contains the communitiy index for
%         node i (where communities are sequentially indexed as 1:M and M
%         is the total number of detected communities).
% 
%         'categories'   ---  a Nx1 array where N is the number of nodes and
%         each entry is either a '0' or '1' denoting the assignment of each 
%         node to one of two communities.

number_nodes=length(partitions);
number_communities=max(partitions);
comm_laterality_array=zeros(number_communities+1,2);
comm_laterality_array(:,1)=[1:number_communities 0];
number_nodes_in_communities=zeros(number_communities,3);

%calcuate the laterality for each community
for i=1:number_communities
    community_i_nodes=find(partitions==i);
    number_nodes_in_i=length(community_i_nodes);
    number_nodes_in_communities(i,1)=number_nodes_in_i;
    categories_i=categories(community_i_nodes);
    nodes_in_comm_0=find(categories_i==0);
    number_nodes_in_comm_0=length(nodes_in_comm_0);
    number_nodes_in_communities(i,2)=number_nodes_in_comm_0;
    nodes_in_comm_1=find(categories_i==1);
    number_nodes_in_comm_1=length(nodes_in_comm_1);
    number_nodes_in_communities(i,3)=number_nodes_in_comm_1;
    laterality_i=abs(number_nodes_in_comm_0-number_nodes_in_comm_1)/number_nodes_in_i;
    comm_laterality_array(i,2)=laterality_i;
end


%calcuate the laterality for the network partition

%need to calcuated "expected" laterality from each community from surrogate
%data
total_nodes_assignments=sum(number_nodes_in_communities,1)
total_number_in_comm_1=total_nodes_assignments(3);
number_surrogates=1000;
surrogate_laterality=zeros(number_communities,number_surrogates);
for j=1:number_surrogates
    randomized_categories=zeros(1,number_nodes);
    rand_perm_assignment=randperm(number_nodes);
    place_in_comm_1=rand_perm_assignment(1:total_number_in_comm_1);
    randomized_categories(place_in_comm_1)=1;
    %now calculate the community laterality for the randomized data
    for i=1:number_communities
        community_i_nodes=find(partitions==i);
        number_nodes_in_i=length(community_i_nodes);
        randomized_categories_i=randomized_categories(community_i_nodes);
        rand_nodes_in_comm_0=find(randomized_categories_i==0);
        rand_number_nodes_in_comm_0=length(rand_nodes_in_comm_0);
        rand_nodes_in_comm_1=find(randomized_categories_i==1);
        rand_number_nodes_in_comm_1=length(rand_nodes_in_comm_1);
        rand_laterality_i=abs(rand_number_nodes_in_comm_0-rand_number_nodes_in_comm_1)/number_nodes_in_i;
        surrogate_laterality(i,j)=rand_laterality_i;
    end
end
expected_comm_laterality=sum(surrogate_laterality,2)/number_surrogates;

network_laterality=(1/number_nodes)*sum(number_nodes_in_communities(:,1).*comm_laterality_array(1:number_communities,2)-expected_comm_laterality);
comm_laterality_array(number_communities+1,2)=network_laterality;

end

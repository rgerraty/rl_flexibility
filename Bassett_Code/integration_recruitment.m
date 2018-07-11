function [Integration, Recruitment] = integration_recruitment(MA, S)
%% Input Module-Allegiance "MA" and community strucutre "S"
%  Output Integration and Recruitment
%%

% transform S to a column vector
if size(S,1) == 1
    S = S'; 
end

num_node = numel(S);
num_cl = max(S);

H = zeros(num_node, max(S));
for i = 1:max(S)
    H(:,i) = S==i;
end
D_H = H'*H;

Recruitment = zeros(num_cl);
Integration = zeros(num_cl);

Recruitment = D_H^-1*H'*MA*H*D_H^-1;
D = diag(diag(Recruitment));
Integration = D^-.5*Recruitment*D^-.5;



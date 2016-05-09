function [i_s,norm_is]=inter_strength(a_mat,labels,g1,g2)

for i=1:size(a_mat,3)
i_s(i)=sum(sum(a_mat(labels==g2,labels==g1,i)))...
		/((sum(labels==g1)*sum(labels==g2)));
norm_is(i)=i_s(i) ...
		/sqrt((sum(sum(a_mat(labels==g2,labels==g2,i)))/((sum(labels==g2)*sum(labels==g2))))...
			*(sum(sum(a_mat(labels==g1,labels==g1,i)))/((sum(labels==g1)*sum(labels==g1)))));

end

end
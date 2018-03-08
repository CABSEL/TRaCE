function acc_w=sort_bounds(gu,gl,z_matrix)

%This function sorts the edges based on upper and lower bounds and z-scores
% from differential gene expression data. 
%Inputs:
% gu= upper bound in matrix form
% gl= lower bound in matrix form
% z_matrix= matrix of z -scores from differential gene expression data
%Output
% acc_w= ranked edges in weighted matrix form


% Last Update 20.05.2014




% Program written by S.M. Minhaz Ud-Dean. (minhazuddean@gmail.com) The
% author accepts no liability for the quality of the information provided
% or for it being correct, complete or up to date. Liability claims against
% the author concerning either material or intellectual damage or other
% detrimental results resulting from the use or non-use of any information
% provided, including any information that is either incomplete or
% incorrect, will therefore be rejected.


acc_w=z_matrix/max(max(z_matrix)); %normalize between zero and one
n=length(gu); % find the number of genes
gT=adj2acc(gu); % calculate the transitive closure of the upper bound
gt=contrex(gl); % calculate the transitive reduction  of the lower bound

add_gT=max(max(((ones(n)-gT)>0).*acc_w));
acc_w=acc_w+add_gT*gT;

add_gu=max(max(((ones(n)-gu)>0).*acc_w)); 
acc_w=acc_w+add_gu*gu;

add_gl=max(max(((ones(n)-gl)>0).*acc_w)); 
acc_w=acc_w+add_gl*gl;

add_gt=max(max(((ones(n)-gt)>0).*acc_w));
acc_w=acc_w+add_gt*gt;

acc_w=acc_w/max(max(acc_w)); %normalize between zero and one


return
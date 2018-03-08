function [z_matrix,acc_z]=analyze_z_score(data_matrix,z_threshold,z_cutoff)
% [z_matrix,acc_z]=analyze_z_score(data_matrix,z_threshold,z_cutoff)
% This function performs z-score transformation of expression data and
% construct the accessibility matrix. See our manuscript titled "Ensemble
% inference and inferability of gene regulatory networks" for detail.
%
% The inputs are as follows:
% data_matrix = a 3D matrix containing expression data (i,j,k) where the
% i-th row corresponds to the KO experiment of gene i and genes in V_KO,
% the j-th column corresponds to gene j, and the k-th aisle corresponds to
% the k-th replicate. See Fig. 2 in the manuscript.
% z_threshold = threshold value for z-score above which a gene is deemed to
% be differentially expressed in a KO experiment.
% z_cutoff = cutoff value (in multiple of standard deviation) for excluding 
% data in the calculation of sample mean and standard deviation.
%
% The outputs are as follows:
% z_matrix= matrix of average z-scores
% acc_z= accessibility matrix

% Last Update 26.05.2014

% Program written by S.M. Minhaz Ud-Dean. (minhazuddean@gmail.com) The
% author accepts no liability for the quality of the information provided
% or for it being correct, complete or up to date. Liability claims against
% the author concerning either material or intellectual damage or other
% detrimental results resulting from the use or non-use of any information
% provided, including any information that is either incomplete or
% incorrect, will therefore be rejected.



num_it=size(data_matrix,1); 
n=size(data_matrix,2);
num_replicate=size(data_matrix,3);
z_matrix=zeros(num_it,n,num_replicate);
for rep=1:num_replicate
    
    data=data_matrix(:,:,rep);
    for i=1:num_it
        data_i=data(:,i);
        
        mu_i=mean(data_i);
        sigma_i=std(data_i);
        
        %%% Exclude outliers %%%
            data_ii=data_i;
            data_ii(abs(data_ii-mu_i)>z_cutoff*sigma_i)=[];
            mu_i=mean(data_ii);
            sigma_i=std(data_ii);
        %%% Exclude outliers %%%
        
        
        for j=1:num_it            
            if sigma_i~=0
                zval=(data_i(j)-mu_i)/sigma_i;
            elseif sigma_i==0
                zval=0;
            end
            z_matrix(j,i,rep)=zval;
        end
    end
end
z_matrix=abs(sum(z_matrix,3))/num_replicate;
z_matrix=remove_diagonal(z_matrix);
acc_z=((z_matrix)>=z_threshold);
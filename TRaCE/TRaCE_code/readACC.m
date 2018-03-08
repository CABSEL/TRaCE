function [ko_acc_matrices,acc_0,ko] = readACC(NUM_ELEM,INPUT_ACC,DESC)
% [ko_acc_matrices,acc_0,ko] = readACC(NUM_ELEM,INPUT_ACC,DESC)
% reads accessibility matrices from a specified file and converts the data
% into input arguments for TRaCE. 
%
% The inputs are as follows:
% NUM_ELEM :  number of genes
% INPUT_ACC : the name of the file containing the accessibility list 
% DESC : the name of the file containing description of the input
% accessibility list
%
% The outputs are as follows:
% ko_acc_matrices: stacked accessibility matrices of G_V_KO(k)
% acc_0: accessibility matrix of the gene regulatory network (GRN) of 
% interest G_null
% ko: list of the genes in V_KO(k) in the same order as the stacking of
% accessibility  matrices in the output ko_acc_matrices
%
% Creation of the input accessibility file (INPUT_ACC): The list contains
% the accessibility data from  G_null and G_V_KO(k) gene regulatory
% networks. Each line shows three numbers (comma separated). The first two
% numbers (i, j) correspond to the index of genes, which should be
% interpreted such that gene j is accessible from gene i. The last number
% should always be written as 1 (following the sparse matrix formatting of
% MATLAB). The segmentation of the list into different accessibility
% matrices is specified in the DESC file. The minimum input is the
% accessibility of the GRN G_null, which must be given at the top of the
% file. Below is an example of the file content:
%
% gene_1,gene_2,1
% ......
% gene_i,gene_j,1
% ........
%
% Creation of the description file (DESC): Each row contains at least two
% numbers (comma separated). The first number corresponds to the line
% number in the INPUT_ACC file, indicating the start of a separate
% accessibility list. The subsequent number(s) corresponds to the indices
% of genes in V_KO(k). The first row is set to (1, 0), indicating
% the accessibility list of the GRN G_null, which is always provided at
% the beginning of the INPUT_ACC file (see above). Below is an example of
% description file:
%
% 1, 0
% Line_2, V_KO1
% .......
% Line_i,V_KOi
% .........

% Last Update 26.05.2014
%
% Program written by S.M. Minhaz Ud-Dean. (minhazuddean@gmail.com). The
% author accepts no liability for the quality of the information provided
% or for it being correct, complete or up to date. Liability claims against
% the author concerning either material or intellectual damage or other
% detrimental results resulting from the use or non-use of any information
% provided, including any information that is either incomplete or
% incorrect, will therefore be rejected.


%% Read and Process Accessibility Matrix of the wild type and KO
ko=DESC(:,2); %separate the list of ko
line_numbers=DESC(:,1); %separate the list of line numbers
num_exp=size(ko,1); %number of accessibility matrices
ko_acc_matrices=zeros(NUM_ELEM,NUM_ELEM, num_exp-1);




for i=1:num_exp
    if i<num_exp
        data_dump=INPUT_ACC(line_numbers(i):line_numbers(i+1)-1,:); 
                               %read required lines from the sparse matrix
    else
        data_dump=INPUT_ACC(line_numbers(i):end,:); 
                    %read required lines from the sparse matrix
    end
    Acc2=spconvert(data_dump); %convert to full matrix
    
    row=size(Acc2,1); %Find the size of the converted full matrix
    column=size(Acc2,2);  %Find the size of the converted full matrix
    Acc1=zeros(NUM_ELEM);
    Acc1(1:row,1:column)=Acc2; 
              %Make a padded accessibility matrix of the appropriate size
    
    lKno=ko(i); %the knocked out gene number
    if (lKno==0)
        acc_0= Acc1;         
    else
        
    ko_acc_matrices(:,:,i-1)=Acc1;
    
    end
end


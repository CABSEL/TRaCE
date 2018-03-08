% An example of the analysis of double knockout expression data using TRaCE
% and z-score analysis

%Inputs
net=1; %Which network of the given examples
z_threshold=2.0; %z_threshold see manuscript for explanation
z_cutoff=3.0; %z_cutoff see manuscript for explanation

% Parameters
num_genes=100; %number of genes
num_ko=100; %number of Single KO experiments
num_replicate=5; %number of replicates

% add the path to trace folder
addpath('../TRaCE_code')

%Declare variables
data_matrix=zeros(num_ko,num_genes,num_replicate);
data_stack=zeros(num_ko,num_genes);

% Read Data Files for Single KO experiments
for rep=1:num_replicate
    file=strcat('./Net',num2str(net),'/replicate',num2str(rep),...
        '/insilico_size100_',num2str(net),'_knockouts.tsv'); %change this 
    % file name and folder according to the location of your data files
    dat=importdata(file,'\t', 1);
    data_matrix(:,:,rep)=dat.data;
    data_stack=data_stack+dat.data;
end

%differential expression analysis for single KO experiments
[z_matrix,acc_z]=analyze_z_score(data_matrix,z_threshold,z_cutoff);

acc_z_tr=contrex(acc_z); 
% bounds 
GU=acc_z; 
GL=acc_z_tr;

% Sort bounds
acc_w=sort_bounds(acc_z,acc_z_tr,z_matrix);
acc_w=acc_w/max(max(acc_w));


% read gold standard file
gold_file=strcat('./Net',num2str(net),'/insilico_size100_',...
    num2str(net),'_goldstandard.tsv');
G=read_gold_standard(num_genes,gold_file);

% display distance of the estimated bounds from the gold standard 
results= [nnz((G-GU)>0)   nnz((GL-G)>0)     nnz((GU~=GL))];
display(results)

% save the ranked edges 
[a,b,s]=find(acc_w);
acc_w_s=[a,b,s];
output_file=strcat('DREAM4_Example_InSilico_Size100_',num2str(net),'.txt');
dlmwrite(output_file,acc_w_s); % This output can be used to calculated AUPR
     %and AUROC


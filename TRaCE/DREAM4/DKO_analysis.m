% An example of the analysis of double knockout expression data using TRaCE
% and z-score analysis

%Inputs
net=1; %Which network of the given examples
z_threshold=2.0; %z_threshold see manuscript for explanation
z_cutoff=3.0; %z_cutoff see manuscript for explanation
threshold=0.65; %threshold for wild-type accessibility matrix estimation 
%see manuscript for explanation


% Parameters
num_genes=100; %number of  genes
num_dko=10000; %number of double KO
num_dko_replicates=5; %number of  replicates for each double KO
num_ko=100; % number of  single KO
num_replicate=5;  %number of replicates for each single KO

% add the path to trace folder
addpath('../TRaCE_code')

%Declare variables
data_matrix=zeros(num_ko,num_genes,num_replicate);
dko_data_matrix=zeros(num_dko,num_genes,num_dko_replicates);
ko_z_matrices=zeros(num_genes,num_genes,num_ko);
ko_acc_z_matrices=zeros(num_genes,num_genes,num_ko);
ko= 0:num_genes; 
% In the examples the double KO data are arranged in the following manner:
% Columns represent the genes
% line 1 represent the KO of gene 1
% line i represent the KO of gene 1 and gene i
% line 101 represent the KO of gene 2 and gene 1
% line 100+i represent the KO of gene 2 and gene i
% line j*100+i represent the KO of gene j+1 and gene i

% Data file locations. Change according to the location of your data files
file_part1=strcat('./Net',num2str(net),'/replicate');
file_part2=strcat('/insilico_size100_',num2str(net),'_knockouts.tsv');
file_part3=strcat('./Net',num2str(net),'/rp');
file_part4=strcat('/insilico_size100_',num2str(net),'_dualknockouts.tsv');

% Read Data Files for Single KO experiments
for rep=1:num_replicate
    file=strcat(file_part1,num2str(rep),file_part2);
    dat=importdata(file,'\t', 1);
    data_matrix(:,:,rep)=dat.data;
    clear dat
end

%differential expression analysis for single KO experiments
[z_matrix,acc_z]=analyze_z_score(data_matrix,z_threshold,z_cutoff);

% Read Data Files for Double KO experiments
for rep=1:num_dko_replicates
    file=strcat(file_part3,num2str(rep),file_part4);
    dat=importdata(file,'\t', 1);
    dko_data_matrix(:,:,rep)=dat.data;
    clear dat
end

%stack differential expression analysis for double KO experiments
for lko=1:num_ko
    data=dko_data_matrix((lko-1)*num_genes+1:lko*num_genes,:,:);
    [z_matrix_l,acc_z_l]=analyze_z_score(data,z_threshold,z_cutoff);
    ko_z_matrices(:,:,lko)=z_matrix_l;
    ko_acc_z_matrices(:,:,lko)=acc_z_l;
   
end



% calculate bounds with TRaCE
[GL,GU]= trace(ko_acc_z_matrices,acc_z,ko,threshold,1);

% sort bounds
acc_w=sort_bounds(GU,GL,z_matrix);
acc_w=acc_w/max(max(acc_w));


% read gold standard file
gold_file=strcat('./Net',num2str(net),'/insilico_size100_',...
    num2str(net),'_goldstandard.tsv');
G=read_gold_standard(num_genes,gold_file);

% display distance of the estimated bounds from the gold standard 
results= [nnz((G-GU)>0)   nnz((GL-G)>0)    nnz((GU~=GL))];
display(results)

% save the ranked edges 
[a,b,s]=find(acc_w);
acc_w_s=[a,b,s];
output_file=strcat('DREAM4_Example_InSilico_Size100_',num2str(net),'.txt');
dlmwrite(output_file,acc_w_s); % This output can be used to calculated AUPR
     %and AUROC


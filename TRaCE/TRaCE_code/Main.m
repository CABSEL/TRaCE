% An application of TRaCE (requires InputNet1.csv and InputIndex1.csv)
%
% Last Update 19.05.2014
%
% Program written by S.M. Minhaz Ud-Dean. (minhazuddean@gmail.com)
% The author accepts no liability for the quality of the information 
% provided or for it being correct, complete or up to date. Liability 
% claims against the author concerning either material or intellectual  
% damage or other detrimental results resulting from the use or non-use of  
% any information provided, including any information that is either 
% incomplete or incorrect, will therefore be rejected.

clear; clc;

%Set the number of genes in the network
num_elem = 100;       

%Load input accessibility matrices
input_acc = dlmread('InputNet1.csv'); 
desc = dlmread('InputIndex1.csv'); 


%Set the parameters of TRACE
threshold=0.75; %threshold in preprocessing step
with_error=0; %input matrices are erroneous

%Define output files
upperbound_file='GU1.csv';
lowerbound_file='GL1.csv';

%Perform TRaCE (with error correction)
[ko_acc_matrices,acc_0,ko]   =readACC(num_elem,input_acc,desc);
[GL,GU]=trace(ko_acc_matrices,acc_0,ko,threshold,with_error);

%Save the output of TRaCE: GU and GL 
GU=double(GU); %Change format to  Write in ASCII
GL=double(GL); % Change format to Write in ASCII
save(upperbound_file, 'GU','-ascii', '-append');
save(lowerbound_file, 'GL','-ascii', '-append');


% Read the true adjacencey
%Load true adjacency matrix for comparison 
true_file='Adj.csv';
adj=dlmread(true_file);
adj=sparse(adj(:,1)',adj(:,2)',adj(:,3)');
adj2=zeros(num_elem);
adj2(1:size(adj,1),1:size(adj,2))=adj;
adj=double(adj2);

%calculate n(GU-G)
GU_G=GU-adj;
GU_G(GU_G<0)=0;
GU_G_d=sum(sum(GU_G));
%calculate n(GL-G)
GL_G=GL-adj;
GL_G(GL_G<0)=0;
GL_G_d=sum(sum(GL_G));
%calculate n(G-GU)
G_GU=adj-GU;
G_GU(G_GU<0)=0;
G_GU_d=sum(sum(G_GU));
%calculate n(G-GL)
G_GL=adj-GL;
G_GL(G_GL<0)=0;
G_GL_d=sum(sum(G_GL));
%Format and display results
results= [ GU_G_d  GL_G_d  G_GU_d  G_GL_d];
display(results)

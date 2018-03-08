function adj=read_gold_standard(n,filename)

% adj=read_gold_standard(n,filename)
% reads the gold standard file from the output file of GeneNetWeaver
% (gnw.sourceforge.net)
%
% The inputs are as follows:
% n: Number of genes
% filename: name of the gold standard file including full path
% 
% The output is as follows:
% adj= adjacency matrix of the gold standard network 

% Last Update 26.05.2014

% Program written by S.M. Minhaz Ud-Dean. (minhazuddean@gmail.com) The
% author accepts no liability for the quality of the information provided
% or for it being correct, complete or up to date. Liability claims against
% the author concerning either material or intellectual damage or other
% detrimental results resulting from the use or non-use of any information
% provided, including any information that is either incomplete or
% incorrect, will therefore be rejected.

fid=fopen(filename); %Open file

adj=zeros(n);
tline = 'll'; %Read First Line
while ischar(tline)  %Continue till end
    tline = fgetl(fid);
    
    x1=tline;
    x1=strrep(x1, '->', ' ');
    
    x1=strrep(x1, 'G', ' ');
    
    
    x1=strrep(x1, '"', ' ');
    try 
    connect = sscanf(x1, '%d');
    catch exception
        break
    end
        
    
    adj(connect(1),connect(2))=connect(3);
    
end

return
function ACC = adj2acc(ADJ)
% ACC = adj2acc(ADJ) returns the transitive closure of the graph defined by
% the adjacency matrix ADJ. The adjacency matrix must be a square matrix.
% 
% Last Update 20.05.2014
%
% Program written by S.M. Minhaz Ud-Dean. (minhazuddean@gmail.com)
% The author accepts no liability for the quality of the information 
% provided or for it being correct, complete or up to date. Liability 
% claims against the author concerning either material or intellectual  
% damage or other detrimental results resulting from the use or non-use of  
% any information provided, including any information that is either 
% incomplete or incorrect, will therefore be rejected.


if (size(ADJ,1)==size(ADJ,2))
    n=length(ADJ);
else
    error('The input must be a square matrix');
end

dummy = expm(double(full(ADJ))) - eye(n);
ACC = (dummy>1e-16);
ACC(ACC<=0)=0;
end
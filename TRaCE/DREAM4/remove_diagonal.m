function b=remove_diagonal(a)
% This function set the diagonal of the matrix a to 0

% Last Update 26.05.2014

for i=1:size(a,1)
    a(i,i)=0;
end
b=a;
end
function ACC_TR=contrex(ACC)
% ACC_TR=contrex(ACC) returns the ConTREx of the input accessibility
% matrix ACC. The accessibility matrix input must be a square matrix.
% 
% Last Update 20.05.2014

% Program written by S.M. Minhaz Ud-Dean. (minhazuddean@gmail.com)
% The author accepts no liability for the quality of the information 
% provided or for it being correct, complete or up to date. Liability 
% claims against the author concerning either material or intellectual  
% damage or other detrimental results resulting from the use or non-use of  
% any information provided, including any information that is either 
% incomplete or incorrect, will therefore be rejected.
%% Condensation of cycles into strong components
num_vertices=size(ACC,1);
mapping=zeros(1,num_vertices); %empty mapping
next_map=0; %No vertex in condensed graph

for i=1:num_vertices
    if (mapping(i)==0) %If i is not already mapped
        next_map=next_map+1; %Get a new vertex in the condensed graph
        mapping(i)=next_map; %Assign that vertex to i
        for j=1:num_vertices
            if (j~=i)%if they are not the same nodes
                if (ACC(i,j)*ACC(j,i)==1)
                    %If i and j are mutually
                    % accessible from each other
                    mapping(j)=next_map; 
                         %Map both on the same element 
                         % in the condensed graph
                end
            end
        end
    end
    
end

mapping2=setdiff(mapping,0); %Remove zeros from the mapping
num_vertices_cond=size(mapping2,2); 
             %Number of unique vertices in the condensed graph
acc_cond=zeros(num_vertices_cond);

%% Adding Appropriate Edges connecting the condensed elements
for i=1:num_vertices
    if (mapping(i)~=0) %Check out empty vertices
        for j=1:num_vertices
            if (j~=i)%if they are not the same nodes
                if(mapping(i)~=mapping(j))  
                    %Check out nodes mapped on the same vertex of 
                    % the condensed graph
                    if(mapping(j)~=0) %Check out empty vertices
                        if (acc_cond(mapping(i),mapping(j))==0) 
                            %If they are not connected in the 
                             % condensed graph
                            if (ACC(i,j)==1) 
                                %But if they are connected in 
                                % the original graph
                                acc_cond(mapping(i),mapping(j))=1;  
                                %Connect them in the condensed graph
                            end
                        end
                    end
                end
            end
        end
    end
end

%% Find the most parsimonious accessibility matrix for
% the condensed accessibility matrix

acc2=double(acc_cond);
acc3=acc2;

for i=1:num_vertices_cond
    acc3=(acc3-boolean(acc2(:,i)*acc2(i,:)));
    acc3(acc3<=0)=0;
end

acc_pars=boolean(acc3);


%% Map it back to the Form of the argument matrix
ACC_TR=zeros(num_vertices); %Initiate with zeros
for i=1:num_vertices_cond %Iterate over the vertices of the condensed graph
    found=find(mapping==i);
    if (size(found,2)==1) %If i is a unique element
        for j=1:num_vertices_cond 
            %Iterate over the vertices of the condensed graph
            if (j~=i)%if they are not the same nodes
                if (size(find(mapping==j),2)==1) %If j is unique
                    if(acc_pars(i,j)==1) 
                        %If they are connected in the condensed 
                         % parsimonious graph
                        ACC_TR(mapping==i,mapping==j)=1; 
                        %Set the corresponding element to 1
                    end
                end
            end
        end
    elseif (size(found,2)==2) %If i is a binary strong element
        j=found(1);%find the positions
        k=found(2);%find the positions
        ACC_TR(j,k)=1; %Set the corresponding element to 1
        ACC_TR(k,j)=1; %Set the corresponding element to 1
    end
end
end
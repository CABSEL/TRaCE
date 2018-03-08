function [GL,GU]=trace(ko_acc_matrices,acc_0,ko,THRESHOLD,ERROR)

% [GL,GU] = trace(ko_acc_matrices,acc_0,ko,THRESHOLD,ERROR)
% Returns the upper and lower bounds.  
%
% The inputs are as follows
% ko_acc_matrices :  stack of KO accessibility matrices
% acc_0 : wild-type (G_O) accessibility matrix
% ko: list of the knockouts in the same order as the accessibility matrices
% in ko_acc_matrices
% THRESHOLD : threshold value for preprocessing step
% ERROR : set to 0 for error-free input, otherwise set to 1.
% 
% The outputs are as follows
% GL: estimated lower bound of the network
% GU: estimated upper bound of the network
%
% Last Update 26.05.2014
%
% Program written by S.M. Minhaz Ud-Dean. (minhazuddean@gmail.com) The
% author accepts no liability for the quality of the information provided
% or for it being correct, complete or up to date. Liability claims against
% the author concerning either material or intellectual damage or other
% detrimental results resulting from the use or non-use of any information
% provided, including any information that is either incomplete or
% incorrect, will therefore be rejected.


NUM_ELEM=size(ko_acc_matrices,1);

num_exp=length(ko);


Kno_Adj_Stack=zeros(NUM_ELEM,NUM_ELEM,num_exp,'int16');
Kno_TR_Stack=zeros(NUM_ELEM,NUM_ELEM,num_exp,'int16');
%Initialize Three Dimensional Array to stack Knockout Accessibility
Kno_TC_Stack=ones(NUM_ELEM,NUM_ELEM,num_exp,'int16'); 
       %consisting of All ones
Kno_TC_Stack_TC=ones(NUM_ELEM,NUM_ELEM,num_exp,'int16'); 
          %consisting of All ones

%Initialize Three Dimensional Array to stack Adjacency Deduced from Knockout
Kno_Adj_Stack2=zeros(NUM_ELEM,NUM_ELEM,num_exp,'int16'); 
                 %consisting of All zeros
Kno_Adj_Stack2_TC=zeros(NUM_ELEM,NUM_ELEM,num_exp,'int16');
                  %consisting of All zeros

for i=1:num_exp

    
    %Use G^t for lower bound
    lKno=ko(i); %the knocked out gene number
    if lKno==0
        Acc1=acc_0;
    else
        Acc1=ko_acc_matrices(:,:,i-1);
    end
    Acc_TC=adj2acc(Acc1); %transitive closure
    Acc_TR=contrex(Acc_TC+0); 
          %transitive reduction of the transitive closure
    Kno_TR_Stack(:,:,i)=Acc_TR; %stack the TR
    
    if (i==1)
        Kno_Adj_Stack(:,:,i)=boolean(int16(Acc_TR));
    else
        Kno_Adj_Stack(:,:,i)=boolean(int16(Acc_TR) ...
            +int16(Kno_Adj_Stack(:,:,i-1)));
    end
    
    
    %Use G^T for upper bound
    Kno_TC_Stack(:,:,i)= Acc1;
    Kno_TC_Stack_TC(:,:,i)= Acc_TC;
    
    if (lKno==0)
        Acc_Wt= Kno_TC_Stack(:,:,i); 
         % The first acc matrix in the input must be the WT(G_O)
        Kno_Adj_Stack2(:,:,i)=Acc_Wt;
        Kno_Adj_Stack2_TC(:,:,i)=Kno_TC_Stack_TC(:,:,i);
    else
        %Matrix to mask (lKno,lKno)
        mask=zeros(NUM_ELEM,NUM_ELEM);
        mask(:,lKno)=1;%Mask row
        mask(lKno,:)=1;%Mask columns
        Kno_Adj_Stack2(:,:,i)=Acc_Wt.*(Kno_TC_Stack(:,:,i)+int16(mask));
        Kno_Adj_Stack2_TC(:,:,i)=Acc_Wt.*(Kno_TC_Stack_TC(:,:,i)...
            +int16(mask));
    end
end

if (ERROR==0)
    %OK
elseif (ERROR==1)
    %OK
else
    error(strcat('ERROR must be set. Set it to 1 if the accessibility',...
        'matrices are erroneous, otherwise set it to 0.'));
end

num_exp=size(Kno_Adj_Stack2,3);
NUM_ELEM=size(Kno_TR_Stack,1);
%%  Error free Reconstruction of Network
if (ERROR==0)
    
    
    Adj_estim=ones(NUM_ELEM,NUM_ELEM,'int16');
    
    
    for l=1:num_exp
        
        
        Adj_estim=Adj_estim.*Kno_Adj_Stack2(:,:,l);
        
    end
    
    GU=double(Adj_estim);
    GL=double(Kno_Adj_Stack(:,:,num_exp));
    return
end

%% =======================================%%
%Error Correction for upper bound from erroneous
%measurements G^m

Acc_Wt_Revised=sum(Kno_TC_Stack,3);
Acc_Wt_Revised(Acc_Wt_Revised<(num_exp*THRESHOLD))=0; 
         %Conversion to binary matrix 
Acc_Wt_Revised(Acc_Wt_Revised>=(num_exp*THRESHOLD))=1; 
       %Conversion to binary matrix 
revised_upperbound=int16(Acc_Wt_Revised);

corrected_lower_bound=int16(contrex(Acc_Wt_Revised)); %initialize with  TR of revised WT

for i_exp=2:num_exp % i_exp=1 for wild-type
    i=ko(i_exp); %find gene number
    mask_testable=int16(Acc_Wt_Revised(:,i)*Acc_Wt_Revised(i,:));
      %Find the edges testable by i-th knockout
    mask_hide_ko=ones(NUM_ELEM,'int16'); 
     %Hiding knockouts. First allow everything to change
    mask_hide_ko(i,:)=0; 
     %Then block the connections from the i-th gene to change
    mask_hide_ko(:,i)=0; 
     % Then block the connections to the i-th gene to change
    mask_allow_change=int16(mask_testable.*mask_hide_ko); 
     %Finally only allow the testable connections neither
     % toward nor from the i-th gene to change
    revised_upperbound=revised_upperbound.*(ones(NUM_ELEM,'int16')...
        -(mask_allow_change-mask_allow_change.*Kno_TC_Stack(:,:,i_exp))); 
             %See justification below
    %justification: mask.*c sets everything other than the unchanged
    %testable edges to zero. Now mask-mask.*c selects the testable
    %edges which changed by this knockout. Then ones(num_elem)-
    %(mask-mask.*c) makes a new mask where everything other than the
    %changed edges are one. Later multiplying it with the bound removes
    %only the changed edges from the bound.
    corrected_lower_bound=int16(boolean(corrected_lower_bound+...
        Kno_TR_Stack(:,:,i_exp).*mask_testable));
end

corrected_upper_bound=revised_upperbound;
corrected_upper_bound(corrected_upper_bound>=1)=1;
   %Conversion to binary matrix 
corrected_upper_bound(corrected_upper_bound<=0)=0;
   %Conversion to binary matrix 


%Find inconsistent edges
% corrected_lower_bound=full(corrected_lower_bound);
errors_found=corrected_lower_bound-corrected_upper_bound;
errors_found(errors_found<0)=0;


%Error Correction for lower bound
for i=1:NUM_ELEM
    for j=1:NUM_ELEM
        if (errors_found(i,j)==1)
            num_against=0;
            num_favor=0;
            for k=1:num_exp
                if (Kno_TR_Stack(i,j,k)==1 && Kno_TC_Stack(i,j,k)==1)
                    num_favor=num_favor+1;
                elseif (Kno_TR_Stack(i,j,k)==0 && Kno_TC_Stack(i,j,k)==0)
                    num_against=num_against+1;
                end
            end
            if num_favor>= num_against
                corrected_upper_bound(i,j)=1;
            elseif num_favor< num_against
                corrected_lower_bound(i,j)=0;
            end
        end
    end
end

%% =======================================%%

%Error Correction for upper bound from TC of erroneous
%measurements TC(G^m)

Acc_Wt_Revised=sum(Kno_TC_Stack_TC,3);
Acc_Wt_Revised(Acc_Wt_Revised<(num_exp*THRESHOLD))=0; 
           %Conversion to binary matrix 
Acc_Wt_Revised(Acc_Wt_Revised>=(num_exp*THRESHOLD))=1; 
                 %Conversion to binary matrix 
revised_upperbound=int16(Acc_Wt_Revised);


for i_exp=2:num_exp  %i_exp=1 for wild-type
    i=ko(i_exp);
    mask_testable=int16(Acc_Wt_Revised(:,i)*Acc_Wt_Revised(i,:)); 
                %Find the edges testable by i-th knockout
    mask_hide_ko=ones(NUM_ELEM,'int16'); 
              %Hiding knockouts. First allow everything to change
    mask_hide_ko(i,:)=0; 
            %Then block the connections from the i-th gene to change
    mask_hide_ko(:,i)=0;
                   % Then block the connections to the i-th gene to change
    mask_allow_change=int16(mask_testable.*mask_hide_ko); 
    %Finally only allow the testable connections neither toward nor
    % from the i-th gene to change
    revised_upperbound=revised_upperbound.*(ones(NUM_ELEM,'int16')-...
        (mask_allow_change-mask_allow_change.*Kno_TC_Stack_TC(:,:,i_exp)));
                            %See justification below
    %justification: mask.*c sets everything other than the unchanged
    %testable edges to zero. Now mask-mask.*c selects the testable
    %edges which changed by this knockout. Then ones(num_elem)-
    %(mask-mask.*c) makes a new mask where everything other than the
    %changed edges are one. Later multiplying it with the bound removes
    %only the changed edges from the bound.
    
end

corrected_upper_bound=revised_upperbound;
corrected_upper_bound(corrected_upper_bound>=1)=1; 
                %Conversion to binary matrix 
corrected_upper_bound(corrected_upper_bound<=0)=0; 
                   %Conversion to binary matrix 



%Find inconsistent edges

errors_found=corrected_lower_bound-corrected_upper_bound;
errors_found(errors_found<0)=0;

%Error Correction for lower bound


for i=1:NUM_ELEM
    for j=1:NUM_ELEM
        if (errors_found(i,j)==1)
            num_against=0;
            num_favor=0;
            for k=1:num_exp
               if (Kno_TR_Stack(i,j,k)==1 && Kno_TC_Stack_TC(i,j,k)==1)
                    num_favor=num_favor+1;
               elseif (Kno_TR_Stack(i,j,k)==0 && Kno_TC_Stack_TC(i,j,k)==0)
                    num_against=num_against+1;
               end
            end
            if num_favor>= num_against
                corrected_upper_bound(i,j)=1;
            elseif num_favor< num_against
                corrected_lower_bound(i,j)=0;
            end
        end
    end
end
if (ERROR==1)
    GL=corrected_lower_bound;
    GU=corrected_upper_bound;
end
GL=double(GL);
GU=double(GU);
end
clear;clc;
load('brainNetSet_SRtest0.02TA.mat'); %Brain network sets based on different net parameters.
nROI=size(brainNetSet{1,9},1); %ROI size nn  
nSubj=length(lab); %# of subject 
p_value=0.01  ; %p-val
 for i=1:nSubj
                originalNet = brainNetSet{1,9}(:,:,i);
				size(originalNet);
                originalNet=originalNet-diag(diag(originalNet)); % remove the non-informative diagal elements
                originalNet=(originalNet+originalNet')/2; % symmetrization
                originalNet=triu(originalNet);
                data(:,i)=reshape(originalNet,nROI^2,1);
            end

 Train_dat = data;

Train_lab = lab;

 
    POSITIVE_data = Train_dat(:,Train_lab==1);
    NEGATIVE_data = Train_dat(:,Train_lab==-1);
    
    % t-test for feature selection
    [tad,p_tmp] = ttest2(double(POSITIVE_data'), double(NEGATIVE_data'));
  % a=find(p_tmp<0.005);
    [a,aa]=sort(p_tmp);
   a(isnan(a))=[];
   num=70;
    a=a(1,1:num);aa=aa(1,1:num);
   % b=p_tmp(a);
   % c=1./b;
    c=1./a;

    e=zeros(1,8100);
    e(aa)=c;
  
   

reshape(e,nROI,nROI);
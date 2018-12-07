clear; clc;
root=cd; addpath(genpath([root '/DATA'])); addpath(genpath([root '/FUN']));
load fMRI80; data=fMRImciNc; clear fMRImciNc;
nSubj=length(lab);
nROI=size(data{1},2);
    ex=-5:5;
    lambda=2.^ex;
  
    nPar=length(lambda);
    brainNetSet=cell(1,nPar);
    m=80*91;n=90;
    k=91;
    ind=0:80:80*91;
    
    q=2;


    
opts=[];

% Starting point
opts.init=2;        % starting from a zero point

% Termination 
opts.tFlag=5;       % run .maxIter iterations
opts.maxIter=100;   % maximum number of iterations

% Normalization
opts.nFlag=0;       % without normalization

% Regularization
opts.rFlag=0;       % the input parameter 'rho' is a ratio in (0, 1)

% Group Property
opts.q=q;           % set the value for q
opts.ind=ind;       % set the group indices

%----------------------- Run the code mtLeastR -----------------------
fprintf('\n mFlag=0, lFlag=0 \n');
opts.mFlag=0;       % treating it as compositive function 
opts.lFlag=0;       % Nemirovski's line search

    for L=1:nPar 
        brainNet=zeros(nROI,nROI,nSubj);
       
        for i=1:nSubj
           
            tmp=data{i};
            tmp=tmp-repmat(mean(tmp')',1,nROI);% centrlization
            data1{i}=tmp;
        end
        sp=zeros(91*80,90);
      
                  for i=1:k
    ind_i=(ind(i)+1):ind(i+1);
   sp(ind_i,:)=data1{i};  
                  end
                  for j=1:nROI
                y=[sp(:,j)];
                A=sp;
                [x, funVal1, ValueL1]= mtLeastR(A, y, lambda(L), opts);
                x(j,:)=0;
                for kk=1:k
                brainNet(:,j,kk) = x(:,kk);
                end
            end
            
        brainNetSet{L}=brainNet;
        fprintf('Done %d/%d networks!\n',L,nPar);
    end
   save('brainNetSet_Groupsparse.mat','brainNetSet','lab');

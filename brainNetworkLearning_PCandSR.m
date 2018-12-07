%% The main procedure for construct brain networks based different methods
% The original ROI time series of all subjects with labels are stored in a mat file.
% For example, the file name of fMRI data is fMRI.mat where each subject
% corresponds to a cell, and in each cell a matrix is used to store the
% time series in its columns one by one. Also, the file include labels and
% subject index.
% method[1] is PC (Pearson's Correlaton);
% method[2] is SR (Sparse Representation based on Jun et al.'s SLEP);
% method[3] is LR (Low-rank Representation based on Jun et al.'s SLEP)
% method[4] is SLR(Sparse and Low-rank Representation rewritten on Jun et al.'s SLEP)
% % Written by Lishan Qiao @UNC, 11/09/2016
%% basic setting, load data, basic statistics
clear; clc;
root=cd; addpath(genpath([root '/DATA'])); addpath(genpath([root '/FUN']));
load fMRI80; data=fMRImciNc; clear fMRImciNc;
nSubj=length(lab);
nROI=size(data{1},2);
method=input('PCC[1],SR[2]:');

%% Network learning based on Pearson's correlation (including sparsification with different thresholds)
if method==1
    lambda=[0 10 20 30 40 50 60 70 80 90 99] % the values lies in [0,100] denoting the sparsity degree
    disp('Press any key:'); pause;
    nPar=length(lambda);
    brainNetSet=cell(1,nPar);
    for L=1:nPar
        brainNet=zeros(nROI,nROI,nSubj);
        for i=1:nSubj
            currentNet=corrcoef(data{i});
            currentNet=currentNet-diag(diag(currentNet));% no link to oneself
            threhold=prctile(abs(currentNet(:)),lambda(L)); % fractile quantile
            currentNet(find(abs(currentNet)<=threhold))=0;
            brainNet(:,:,i)=currentNet;
        end
        brainNetSet{L}=brainNet;
        fprintf('Done %d/%d networks!\n',L,nPar);
    end
    save('brainNetSet_PC.mat','brainNetSet','lab');
end

%% Network learning based on sparse representation(SR) - SLEP
if method==2
    %Parameter setting for SLEP
    ex=-5:5;
    lambda=2.^ex
    disp('Press any key:'); pause;
    nPar=length(lambda);
    brainNetSet=cell(1,nPar);
    
    opts=[];
    opts.init=2;% Starting point: starting from a zero point here
    opts.tFlag=0;% termination criterion
    % abs( funVal(i)- funVal(i-1) ) ¡Ü .tol=10e?4 (default)
    %For the tFlag parameter which has 6 different termination criterion.
    % 0 ? abs( funVal(i)- funVal(i-1) ) ¡Ü .tol.
    % 1 ? abs( funVal(i)- funVal(i-1) ) ¡Ü .tol ¡Á max(funVal(i),1).
    % 2 ? funVal(i) ¡Ü .tol.
    % 3 ? kxi ? xi?1k2 ¡Ü .tol.
    % 4 ? kxi ? xi?1k2 ¡Ü .tol ¡Á max(||xi||_2, 1).
    % 5 ? Run the code for .maxIter iterations.
    opts.nFlag=0;% normalization option: 0-without normalization
    opts.rFlag=0;% regularization % the input parameter 'rho' is a ratio in (0, 1)
    opts.rsL2=0; % the squared two norm term in min  1/2 || A x - y||^2 + 1/2 rsL2 * ||x||_2^2 + z * ||x||_1
    fprintf('\n mFlag=0, lFlag=0 \n');
    opts.mFlag=0;% treating it as compositive function
    opts.lFlag=0;% Nemirovski's line search
    
    for L=1:nPar
        brainNet=zeros(nROI,nROI,nSubj);
        for i=1:nSubj
            tmp=data{i};
            tmp=tmp-repmat(mean(tmp')',1,nROI);% data centralization
            currentNet=zeros(nROI,nROI);
            for j=1:nROI
                y=[tmp(:,j)];
                A=[tmp(:,setdiff(1:nROI,j))];
                [x, funVal1, ValueL1]= LeastR(A, y, lambda(L), opts);
                currentNet(setdiff(1:nROI,j),j) = x;
            end
            brainNet(:,:,i)=currentNet;
        end
        brainNetSet{L}=brainNet;
        fprintf('Done %d/%d networks!\n',L,nPar);
    end
    save('brainNetSet_SR.mat','brainNetSet','lab');
end



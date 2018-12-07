clear
close all

load brainNetSet_SRtest0.02TA.mat; %Brain network sets based on different net parameters.
root=cd; addpath(genpath([root '/DATA'])); addpath(genpath([root '/FUN']));

nROI=size(brainNetSet{1,1},1); %ROI size
nSubj=length(lab); %# of subject
[nL nZ]=size(brainNetSet); % # of candidate parameter Lambda1 and Lambda2
p_value=0.01; %p-value of t-test for feature selection

kfoldin=nSubj; % inner LOO CV



 

accnum=0;tp=0;tn=0;

Vali_res=zeros(nL,nZ); % for results on validation set (split from training set)
 for iL=1:nL % corresponds to regularized parameter lambda 2 in our model
   for iZ=1:nZ % lambda 1 (L1-norm regularized parameter in our model)

            data=zeros(nROI^2,nSubj);
        
            for i=1:nSubj
                originalNet = brainNetSet{iL,iZ}(:,:,i);
                originalNet=originalNet-diag(diag(originalNet)); % remove the non-informative diagal elements
                originalNet=(originalNet+originalNet')/2; % symmetrization
                originalNet=triu(originalNet);
                data(:,i)=reshape(originalNet,nROI^2,1);
            end
            idx=find(sum(abs(data'))>1e-12); % for removing the trivial zero-rows
            data=data(idx,:);
            
         
            
         
            c_in = cvpartition(91,'k',kfoldin); % leave one out for validation (i.e., parameter slection)
            acc=0; % for save current validation result
            for fdin=1:c_in.NumTestSets
           
                % partion the training samples into a training subset and a
                % validation subset for selecing optimal regularized network parameters.
                InTrain_dat = data(:,training(c_in,fdin));
                InTrain_lab = lab(training(c_in,fdin));
                Vali_dat = data(:,test(c_in,fdin));
                Vali_lab = lab(test(c_in,fdin));
                
                POSITIVE_data = InTrain_dat(:,InTrain_lab==1);
                NEGATIVE_data = InTrain_dat(:,InTrain_lab==-1);
                
                % t-test for feature selection
                [tad,p_tmp] = ttest2(double(POSITIVE_data'), double(NEGATIVE_data'));
                [ordered_p,indp]=sort(p_tmp);
                index = indp(ordered_p<p_value);
                if size(index,2)==0
                    fprintf('invalid=%d,from%d\n',fdin,iZ);
                end
                
                InTrain_dat = InTrain_dat(index,:);
                Vali_dat = Vali_dat(index,:);

                % scaling the data, which is recommended by LIBSVM toolbox for classification.
                MaxV=(max(InTrain_dat'))';
                MinV=(min(InTrain_dat'))';
                [R,C]= size(InTrain_dat);
                InTrain_dat=2*(InTrain_dat-repmat(MinV,1,C))./(repmat(MaxV,1,C)-repmat(MinV,1,C))-1;
                [R,C]= size(Vali_dat);
                Vali_dat=2*(Vali_dat-repmat(MinV,1,C))./(repmat(MaxV,1,C)-repmat(MinV,1,C))-1;
                
                cmd = ['-t 0 -c 1 -q'];% linear kernel
                model = svmtrain(InTrain_lab', InTrain_dat', cmd); % Linear Kernel
                
                [predict_label, accuracy_svm,prob_estimates ] = svmpredict(Vali_lab', Vali_dat', model, '-q');
                
                acc = acc+sum(predict_label==Vali_lab')/length(Vali_lab);
               % fprintf('current foldin=%d, Z=%d, Lambda=%d, \n',fdin,iZ,iL);
            end
            Vali_res(iL,iZ)=acc/c_in.NumTestSets;
       end
    end

    [row col]=find(Vali_res==max(Vali_res(:)));% find optimal regularized network parameters
    lMax=row(end); zMax=col(end);
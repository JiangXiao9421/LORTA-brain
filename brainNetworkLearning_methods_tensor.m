addpath(genpath(cd))
clear


load brainNetSet_SR.mat; 
for j=1:11


X = brainNetSet{j};
  

%maxP = max(abs(X(:)));
%[n1,n2,n3] = size(X);
%Xn = X;
%rhos = 0.3
%ind = find(rand(n1*n2*n3,1)<rhos);
%Xn(ind) = rand(length(ind),1);

 opts = [];

%% RPCA
% lambda = 1/(sqrt(max(n1,n2)));
% Xhat = zeros(n1,n2,n3);
% Shat = Xhat;
% 
% tol = 1e-6;
% maxIter = 500;
% mu = 1e-3;
% for j = 1 : 3
%     Xni = Xn(:,:,j);
%     [Xhat(:,:,j),S,obj,err,iter] = rpca(Xni,lambda,opts);
% end
% Xhat = max(Xhat,0);
% Xhat = min(Xhat,maxP);
% Shat = max(Shat,0);
% Shat = min(Shat,maxP);
% Lr_RPCA = norm(X(:)-Xhat(:))/norm(X(:));
% psnr_RPCA = PSNR(X,Xhat,maxP)
% 
% 
% figure(2)
% subplot(1,3,1)
% imshow(X/max(X(:)))
% subplot(1,3,2)
% imshow(Xn/max(Xn(:)))
% subplot(1,3,3)
% imshow(Xhat/max(Xhat(:)))

%% Tensor RRPCA based on SNN
% alpha = [15 15 1.5];
%  
% [Xhat,E,err,iter] = trpca_snn(Xn,alpha,opts);
% 
% err
% iter
%  
% Xhat = max(Xhat,0);
% Xhat = min(Xhat,maxP);
% psnr = PSNR(X,Xhat,maxP)
% 
% figure(1)
% subplot(1,3,1)
% imshow(X/max(X(:)))
% subplot(1,3,2)
% imshow(Xn/max(Xn(:)))
% subplot(1,3,3)
% imshow(Xhat/max(Xhat(:)))
% 
% pause 

%% Tensor RRPCA based on TNN
[n1,n2,n3] = size(X);
%lambda = 1/sqrt(max(n1,n2)*n3);
lambda = 0.02;
[Xhat,E,err,iter] = trpca_tnn(X,lambda,opts);
for i=1:91
    Xhat(:,:,i)=Xhat(:,:,i)-diag(diag(Xhat(:,:,i)));
end
brainNetSet{j}=Xhat;

 fprintf('Done %d/%d networks!\n',j,11);
%end
 save('brainNetSet_SRtest0.02TA.mat','brainNetSet','lab');





clear all
close all

addpath(genpath('../UGM/'))
addpath(genpath('../MRF/'))
addpath(genpath('../hybrid/'))

%% Algorithm constants
lambda = 3;
beta = 1.5;
Epsilon = 3;
superpixel = 1;

% Load prior
load('../data/prior/embPrior.mat');

%% Prior Pre-processing
Phi = 1 - Phi;
Phi(Phi>0.5) = 1;
Phi(Phi<=0.5) =0;
prior_phis = Phi(:,:,:);
no_ex = size(prior_phis,3);
% Create sd_phis
for k = 1:no_ex
    sd_phis(:,:,k) = make_sdfunc(prior_phis(:,:,k));
end

% Compute kernel sigma
%kernel_sigma = mean(get_Hdistance(circshift(prior_phis,[0 0 -1]),prior_phis,Epsilon));
% kernel sigma or kernel sigma square?
Sigma = zeros(no_ex,no_ex);
whichHeaviside = 'tan';
for i = 1:(no_ex-1)
  for j = (i+1):no_ex
    diff_phi_temp = Heaviside(sd_phis(:,:,i), Epsilon,whichHeaviside) - Heaviside(sd_phis(:,:,j), Epsilon,whichHeaviside);
    Sigma(i,j) = sum(sum(diff_phi_temp.^2));
    Sigma(j,i) = Sigma(i,j);
  end
end

%figure; imagesc(Sigma); colorbar;
kernel_sigma = sum(Sigma(:))/(no_ex^2 - no_ex);

%kernel_sigma=1e-2;
no_ex = size(prior_phis,3);

% loadData
%I = imread('../data/embryo.png');
superpixel = 1; loadDataFruitfly; I = im2double(X);
[nRows,nCols,dump] = size(I);

% Make priors and images match for CV
prior_size = size(Phi(:,:,1));
I_CV = imresize(mean(I,3),prior_size);
I_CV = 1 - I_CV;
I_CV = I_CV/max(I_CV(:));

%% Superpixels
%{
rng(215)
createSuperPixels;
I = 1-mean(I,3);
I = I/max(I(:));
nSp = length(unique(Sp(:)));
%}

doSuperPixels = false;
if doSuperPixels
  rng(215)
  createSuperPixels;
  numNbThreshold = 10;
  nStates = 2;
  createNeighbor;
  mergeSmallSuperPixels;
  numNbThreshold = 10;
  createNeighbor;    
  nNodes = nSp;
end
load([imgName,'sp.mat']);

% Resize superpixels to ChanVese
Sp_r = imresize(Sp,prior_size,'nearest');
nSp = length(unique(Sp_r));
idx0 = unique(Sp_r);
for i = 1:nSp
    Sp_r(find(Sp_r == idx0(i))) = i;
end
Sp = imresize(Sp_r,size(I(:,:,1)),'nearest');
nStates = 2;

numNbThreshold = 10; 
createNeighbor;
[edgePot,edgePotVec] = produceEdgePotentials(par.pb_emag,edgeStruct,idxArray,nStates);
plotSuperPixelNeighbors;

pause
G = zeros(nSp,2);  % index matrix for centers of gravity in CV image size
no_sp = zeros(nSp,1);
for i = 1:nSp
    [row_i col_i] = find(Sp_r == i);
    ind{i} = find(Sp_r == i);
    no_sp(i) = length(ind{i});
    G(i,:) = ceil(mean([row_i col_i],1));
end
%figure; imagesc(Sp_r);
%text(G(:,2), G(:,1), num2strBatch(1:nSp));
idxTemp = G(:,1)+(G(:,2)-1)*prior_size(1);

% transform the pixel intensity to superpixel intensity

X = zeros(nSp,2);
for i = 1:nSp
  X(i,1) = 1-mean(I(find(Sp == i)));
  X(i,2) = std(I(find(Sp == i)));
end
X(:,1) = X(:,1)/max(X(:,1));

%{
figure;
imagesc(labelSuperpixelsOnImage(Sp,X(:,1))); colorbar;

figure;
imagesc(labelSuperpixelsOnImage(Sp,X(:,2))); colorbar;

pause
%}

%% Level set function initialization
m = zeros(prior_size);
m(10:size(I_CV,1)-10,5:size(I_CV,2)-10) = 1;
%  m0 = Phi(:,:,23);
%  m(m0>0.5) = 1;
%  m(m0<=0.5) = 0;
phi = make_sdfunc(m); %m; %(:); %

%% Algorithm
% MRF Initialization
%MuEst = [0.4,.5]; %0.41, 0.75
%SigmaEst = [0.13,0.119];
if size(X,2) ==2
  MuEst = [0.6,0.7; 0.01,0.05];
  SigmaEst = [0.1,0.1; 0.0045,0.0092];

  %MuEst = [0.03,0.76; 0.03,0.1];
  %SigmaEst = [0.032,0.2; 0.03,0.018];


else 
  MuEst = [0,.3];
  SigmaEst = [.1,.1];
end

i = 0;
diff = 1;
iter_max = 20;
nodeBelMF_CV = zeros(prior_size(1)*prior_size(2),nStates);
probLogistic_CV =zeros(prior_size(1)*prior_size(2),1);

size_ICV = size(I_CV);
phi_vec = phi(:);
expLambdaPhi = exp(lambda*phi(idxTemp));
probLogistic = 1./(expLambdaPhi + 1);
probA = ones(nSp,nStates)*(1/nStates);
%probA = [1-probLogistic,probLogistic];%
%probA = [probLogistic,1-probLogistic];
mask = phi;
for i = 1:10
    
    
    %probA = ones(nSp,nStates)*(1/nStates);
    %spIdx = [5,25,24,50,23,49];
    spIdx = [14,41,33,1];
    for s = spIdx	
 %     probA(s,:) = [1,0];
    end
    
    figure;
    imagesc(labelSuperpixelsOnImage(Sp_r,probA(:,2)));
    title('prob(y_i=1|a_i)');
    colorbar
    
    
    for j = 1:10
    %[MuEst,SigmaEst] = normalMixtureEM(MuEst,SigmaEst,X,edgeStruct,probA,20);
    %[nodePot,edgePot] = producePotentials(X,edgeStruct,MuEst,SigmaEst,probA);
      nodePot = producePotentials(X,edgeStruct,MuEst,SigmaEst,probA);        
      [nodeBelMF,edgeBelMF,logZMF] = UGM_Infer_MeanField(nodePot,edgePot,edgeStruct);
      %[nodeBelMF,edgeBelMF,logZMF] = UGM_Infer_LBP(nodePot,edgePot,edgeStruct);    
      [MuEst,SigmaEst] = em_max(nodeBelMF,X);
    
    end

    figure;
    imagesc(labelSuperpixelsOnImage(Sp,nodePot(:,2))); colorbar;
    title('Node potential (y=1)');
    

    figure;
    imagesc(labelSuperpixelsOnImage(Sp,nodeBelMF(:,2))); colorbar;
    title('Mean Field Estimates of Marginals (y=1)');
    pause
    
    % Convert into CV
    for j = 1:nSp;
        nodeBelMF_CV(ind{j},:) = repmat(nodeBelMF(j,:),[no_sp(j) 1]);
        %probLogistic_CV(ind{j}) = probLogistic(j);
    end
%     EQ_M = lambda*(probLogistic_CV.*nodeBelMF_CV(:,2) - (1-probLogistic_CV).*nodeBelMF_CV(:,1));
%     EQ_M = reshape(EQ_M,size_ICV(1),size_ICV(2));
    %m = chanvese_exp(I_CV,m,150,0.1,beta*EQ_M,1,prior_phis, kernel_sigma, Epsilon);
    
    % Take the most likely binary mask of MRF as initialization
    
    %if i>=2  
    %  mask = getmask_frompot(nodeBelMF_CV,size_ICV(1),size_ICV(2));
    %  mask = make_sdfunc(mask);
    %end
    mask = phi;
    %kernel_sigma = 0;
    hybrid_b = 1;
    mu = .1; % smoothness

    
    %kernel_sigma = 1;
    m = chanvese_shiftInvar(I_CV,mask,500,mu,nodeBelMF_CV,hybrid_b,prior_phis(:,:,:), kernel_sigma, Epsilon,beta,1);
    
				%ImgTemp = imgShift(prior_phis(:,:,1),0,20);
				%ImgTemp(:,80:82) = 1;
				%ImgTemp(1:10,:) = 1;
    
				%ImgTemp = ImgTemp + 0.3*randn(size(ImgTemp));
				%    figure;imagesc(ImgTemp);
    
				%m = chanvese_shiftInvar(ImgTemp,imgShift(prior_phis(:,:,2),20,10),1000,mu,nodeBelMF_CV,hybrid_b,prior_phis(:,:,1), kernel_sigma, Epsilon,beta,1);
    
    phi = make_sdfunc(m);   % phi is a matrix!
    % Convert into MRF Superpixel stuff
    phi_vec = phi(:);
    expLambdaPhi = exp(lambda*phi(idxTemp));
    probLogistic = 1./(expLambdaPhi + 1);
    %probA = [1-probLogistic, probLogistic];
    probA = [probLogistic, 1-probLogistic];
    %pause
    
    bdd = findboundary(m,2);
    figure;imshow(bdd+I_CV);
    
end


bdd2 = imresize(bdd,[size(I,1),size(I,2)]);
figure;imshow(I(:,:,1)-bdd2);



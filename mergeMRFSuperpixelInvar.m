clear all
close all

addpath(genpath('../UGM/'))
addpath(genpath('../MRF/'))
addpath(genpath('../Chan-vese'))
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


% loadData
%I = imread('../data/embryo.png');
superpixel = 1; loadDataFruitfly; I = im2double(X);
[nRows,nCols,dump] = size(I);

% Make priors and images match for CV
prior_size = size(Phi(:,:,1));
I_CV = imresize(mean(I,3),prior_size);
I_CV = 1 - I_CV;
I_CV = I_CV/max(I_CV(:));


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
idxSP = G(:,1)+(G(:,2)-1)*prior_size(1);

% transform the pixel intensity to superpixel intensity

X = zeros(nSp,2);
for i = 1:nSp
  X(i,1) = 1-mean(I(find(Sp == i)));
  X(i,2) = std(I(find(Sp == i)));
end
X(:,1) = X(:,1)/max(X(:,1));


figure;
imagesc(labelSuperpixelsOnImage(Sp,X(:,1))); colorbar;

figure;
imagesc(labelSuperpixelsOnImage(Sp,X(:,2))); colorbar;

pause


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
expLambdaPhi = exp(lambda*phi(idxSP));
probLogistic = 1./(expLambdaPhi + 1);
probA = ones(nSp,nStates)*(1/nStates);
%probA = [1-probLogistic,probLogistic];%
%probA = [probLogistic,1-probLogistic];
mask = phi;
for i = 1:10
    
    figure;
    imagesc(labelSuperpixelsOnImage(Sp_r,probA(:,2)));
    title('prob(y_i=1|a_i)');
    colorbar
    
    
    for j = 1:10
      nodePot = producePotentials(X,edgeStruct,MuEst,SigmaEst,probA);        
      [nodeBelMF,edgeBelMF,logZMF] = UGM_Infer_MeanField(nodePot,edgePot,edgeStruct);
      [MuEst,SigmaEst] = em_max(nodeBelMF,X);
    end

    figure;
    temp = log(nodePot(:,2)./nodePot(:,1));    
    imagesc(labelSuperpixelsOnImage(Sp,temp)); colorbar;
    title('Node potential (y=1)');    

    figure;
    imagesc(labelSuperpixelsOnImage(Sp,nodeBelMF(:,2))); colorbar;
    title('Mean Field Estimates of Marginals (y=1)');
    pause
    
    % Convert into CV
    for j = 1:nSp
        nodeBelMF_CV(ind{j},:) = repmat(nodeBelMF(j,:),[no_sp(j) 1]);
    end
    
    mask = phi;
    mu = 1; % smoothness
    Epsilon = 3;
    beta = 1.5;
    lambda1 = 1;
    lambda2 = 3;
    priorWeight = 1;
    if i == 1
      [m,phi_crap,param] = chanvese_invar(I_CV,mask,200,mu, nodeBelMF_CV, prior_phis(:,:,1:3), kernel_sigma, Epsilon, beta, lambda1,priorWeight);
    else
      [m,phi_crap,param] = chanvese_invar(I_CV,mask,200,mu, nodeBelMF_CV, prior_phis(:,:,1:3), kernel_sigma, Epsilon, beta, lambda1,priorWeight,param);
    end
    
    
    phi = make_sdfunc(m);   % phi is a matrix!
    % Convert into MRF Superpixel stuff
    %phi_vec = phi(:);
    expLambdaPhi = exp(lambda2*phi(idxSP));
    probLogistic = 1./(expLambdaPhi + 1);
    probA = [probLogistic, 1-probLogistic];
    
    %figure;
    %imagesc(labelSuperpixelsOnImage(Sp_r,probA(:,2)));
   
    bdd = findboundary(m,2);
    figure;imshow(bdd+I_CV);
    
end


bdd2 = imresize(bdd,[size(I,1),size(I,2)]);
figure;imshow(I(:,:,1)-bdd2);



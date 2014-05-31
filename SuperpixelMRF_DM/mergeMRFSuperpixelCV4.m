clear all
close all

tic
addpath(genpath('../UGM/'))
addpath(genpath('../MRF'))

%% Algorithm constants

% Choosing which combination to run
% kernel_sigma = 0 -> no priors
% hybrid_b = 0 -> only Level set
priors_b = 1;
hybrid_b = 1;

% Choosing some constants
lambda = 3; % 5, beta 1.5%
beta = 1.5;
Epsilon = 3;

% Run what
superpixel = 1;
plot_all = 0; % 0 if you want results only
plot_int_results = 0; % 1 if you want intermediate results only
doSuperPixels = 0; %  if new data, set this 1 to run superpixel algorithm

whichHeaviside = 'sin';

%% Prior Pre-processing
% Load prior
load('../data/prior/shapePriorShift.mat');

Phi = 1 - Phi;
Phi(Phi>0.5) = 1;
Phi(Phi<=0.5) =0;
prior_phis = Phi(:,:,:);
no_ex = size(prior_phis,3);

% Compute kernel sigma
if priors_b == 1
    kernel_sigma = mean(get_Hdistance(circshift(prior_phis,[0 0 -1]),prior_phis,Epsilon,whichHeaviside));
else
    kernel_sigma=0;
end

% Create sd_phis
for k = 1:no_ex
    sd_phis(:,:,k) = make_sdfunc(prior_phis(:,:,k));
end


%% Load Data
loadDataFruitfly; I = im2double(X);
[nRows,nCols,dump] = size(I);

% Make priors and images match for CV
prior_size = size(Phi(:,:,1));
I_CV = imresize(mean(I,3),prior_size);
I_CV = 1 - I_CV;
I_CV = I_CV/max(I_CV(:));

%% Superpixels

% Run superpixel lgorithm
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
    save([imgName,'sp.mat'],'Sp','nSp','nNodes','edgeStruct','idxArray','nStates','nRows','nCols','par','I_sp','numNbThreshold');
end
load([imgName,'sp.mat']);
toc

% Resize superpixels to ChanVese
Sp_r = imresize(Sp,prior_size,'nearest');
nSp = length(unique(Sp_r));
idx0 = unique(Sp_r);
for i = 1:nSp
    Sp_r(find(Sp_r == idx0(i))) = i;
end
Sp = imresize(Sp_r,size(I(:,:,1)),'nearest');
nStates = 2;

% Create edge structure
createNeighbor;
[edgePot,edgePotVec] = produceEdgePotentials(par.pb_emag,edgeStruct,idxArray,nStates);
plotSuperPixelNeighbors;
%pause

% Create index matrix for centers of gravity in CV image size
G = zeros(nSp,2);  
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

% Transform the pixel intensity to superpixel intensity
X = zeros(nSp,2);
for i = 1:nSp
    X(i,1) = 1-mean(I(find(Sp == i)));
    X(i,2) = std(I(find(Sp == i)));
end
X(:,1) = X(:,1)/max(X(:,1));

if plot_all ==1
    figure;
    imagesc(labelSuperpixelsOnImage(Sp,X(:,1))); colorbar;
    
    figure;
    imagesc(labelSuperpixelsOnImage(Sp,X(:,2))); colorbar;
    
    pause
end

%% Level set function initialization
m = zeros(prior_size);
m(10:size(I_CV,1)-10,5:size(I_CV,2)-10) = 1;
% % --- Use different initialization
%  m0 = Phi(:,:,23);
%  m(m0>0.5) = 1;
%  m(m0<=0.5) = 0;
phi = make_sdfunc(m);

%% Algorithm

% MRF Initialization
%MuEst = [0.4,.5]; %0.41, 0.75
%SigmaEst = [0.13,0.119];
if size(X,2) ==2
    MuEst = [0.6,0.7; 0.01,0.05];
    SigmaEst = [0.1,0.1; 0.0045,0.0092];
else
    MuEst = [0,.3];
    SigmaEst = [.1,.1];
end

i = 0;
diff = 1;
iter_max = 20;

% Initialize potentials and contour
nodeBelMF_CV = zeros(prior_size(1)*prior_size(2),nStates);
probLogistic_CV = zeros(prior_size(1)*prior_size(2),1);
size_ICV = size(I_CV);
phi_vec = phi(:);
probA = ones(nSp,nStates)*(1/nStates);

% % Initializing with initial contour
% expLambdaPhi = exp(lambda*phi(idxTemp));
% probLogistic = 1./(expLambdaPhi + 1);
%probA = [probLogistic,1-probLogistic];     

mask = phi;


for i = 1:4
    
    %% Run EM for parameter estimation & Q_M^i calculation
    for j = 1:10
        nodePot = producePotentials(X,edgeStruct,MuEst,SigmaEst,probA);
        [nodeBelMF,edgeBelMF,logZMF] = UGM_Infer_MeanField(nodePot,edgePot,edgeStruct);
        [MuEst,SigmaEst] = em_max(nodeBelMF,X);
    end
    
    
    % Convert into CV
    for j = 1:nSp;
        nodeBelMF_CV(ind{j},:) = repmat(nodeBelMF(j,:),[no_sp(j) 1]);
        %probLogistic_CV(ind{j}) = probLogistic(j);
    end
    
    
    %% Plotting for testing reasons
    
    if plot_all == 1
        figure;
        imagesc(labelSuperpixelsOnImage(Sp_r,probA(:,2)));
        title('prob(y_i=1|a_i)');
        colorbar
        
        figure;
        imagesc(labelSuperpixelsOnImage(Sp,nodePot(:,2))); colorbar;
        title('Node potential (y=1)');
        
        figure;
        imagesc(labelSuperpixelsOnImage(Sp,nodeBelMF(:,2))); colorbar;
        title('Mean Field Estimates of Marginals (y=1)');
        pause
    end
    
    %% Level set segmentation for calculation of p(y_i|a_i)
    % % Take the most likely binary mask of MRF as initialization
    %if i>=2
    %  mask = getmask_frompot(nodeBelMF_CV,size_ICV(1),size_ICV(2));
    %  mask = make_sdfunc(mask);
    %end
    
    % Run Level Set
    mask = phi;
    [phi_crap m] = chanvese_exp(I_CV,mask,100,0.1,nodeBelMF_CV,hybrid_b,sd_phis, kernel_sigma, Epsilon,beta,1,whichHeaviside);
    % Make it a signed distance func
    phi = make_sdfunc(m);   % phi is a matrix!
    
    % Compute p(y_i|a_i)
    phi_vec = phi(:);
    expLambdaPhi = exp(lambda*phi(idxTemp));
    probLogistic = 1./(expLambdaPhi + 1);
    probA = [probLogistic, 1-probLogistic];
    
    %% Plotting intermediate results
    bdd = findboundary(m,2); %figure;imshow(bdd+I_CV);
    if plot_int_results == 1 || plot_all ==1
        bdd2 = imresize(bdd,[size(I,1),size(I,2)]);
        figure;imshow(I(:,:,1)-bdd2);
        
        % LevelSet result
        help_mat = reshape(nodeBelMF_CV(:,2),[size(I_CV,1),size(I_CV,2)]);
        help_mat = get_binarymask(help_mat);
        bdd_mrf = findboundary(help_mat,2);
        bdd_mrf = imresize(bdd_mrf,[size(I,1),size(I,2)]);
        figure; imshow(I(:,:,1) - bdd_mrf);
    end
    
end

%% Plot final results

% ChanVese result
bdd2 = imresize(bdd,[size(I,1),size(I,2)]);
figure;imshow(I(:,:,1)-bdd2);

% LevelSet result
help_mat = reshape(nodeBelMF_CV(:,2),[size(I_CV,1),size(I_CV,2)]);
help_mat = get_binarymask(help_mat);
bdd_mrf = findboundary(help_mat,2);
bdd_mrf = imresize(bdd_mrf,[size(I,1),size(I,2)]);
figure; imshow(I(:,:,1) - bdd_mrf);


toc

clear all
close all

addpath('../MRF/');
% addpath('../Chan-Vese/');


%% Algorithm constants
lambda = 1;
beta = 1;
Epsilon = 1;
superpixel = 0;


%% Load Data
load('../data/prior/shapePriorShift.mat');loadDataFruitfly;
%generateData;
% Create graph
edgeStruct = createImgStructure(X,nStates);



%% Prior Pre-processing
Phi = 1-Phi;
Phi(Phi>0.5) = 1;
Phi(Phi<=0.5) =0;
prior_phis = Phi(:,:,[1,4,9,11,16,20,26]);
%prior_phis = Phi(:,:,[1,20]); 
no_ex = size(prior_phis,3);

% For speeding up:
% Calculate the diffmat between the first and the rest. 
% Then only need to compute diffmat betwen current and the first
%diff_phi_mat = Heaviside(repmat(prior_phis(:,:,1), [1 1 no_ex-1]), Epsilon) - Heaviside(sd_phis, Epsilon);


% Compute kernel sigma
kernel_sigma = mean(get_Hdistance(circshift(prior_phis,[0 0 -1]),prior_phis,Epsilon));
%kernel_sigma = 1;

%% Level set function initialization
m = zeros(size(X,1),size(X,2));
m(10:size(X,1)-10,5:size(X,2)-10) = 1;
% m0 = Phi(:,:,36);
% m(m0>0.5) = 1;
% m(m0<=0.5) = 0;
phi = m(:); %make_sdfunc(m);

% %% MRF Intialization
% MuEst = [0,1];
% SigmaEst = [0.1,0.1];
%   
% for i = 1:10
%   expLambdaPhi = exp(lambda*phi(:));
%   probLogistic = 1./(expLambdaPhi + 1);
%   probA = [probLogistic, 1 - probLogistic];
%   
% %   figure;
% %   imagesc(reshape(probA(:,2)./probA(:,1),nRows,nCols));
% %   title('prob(y_i|a_i)');
% %   colorbar
% %   
%   %pause
% 
%   %nodePot_old = nodePot;
%   [MuEst,SigmaEst] = normalMixtureEM(MuEst,SigmaEst,X,edgeStruct,probA,5);
%   [nodePot,edgePot] = producePotentials(X,edgeStruct,MuEst,SigmaEst,probA);
%   nodePot_norm = sum(nodePot,2);
%   
% %   figure;
% %   imagesc(reshape(nodePot(:,2)./nodePot_norm,nRows,nCols));
% %   title('Node Potential Difference');  
% %   colorbar
% %   pause
%   
%   [nodeBelMF,edgeBelMF,logZMF] = UGM_Infer_MeanField(nodePot,edgePot,edgeStruct);  
%   
%   figure;
%   imagesc(reshape(nodeBelMF(:,2),nRows,nCols));
%   title('Mean Field Estimates of Marginals');  
%   colorbar
%   %pause
% 
% %     EQ_M = lambda*(probLogistic.*nodeBelMF(:,2) - (1-probLogistic).*nodeBelMF(:,1));
% %     EQ_M = reshape(EQ_M,size(X,1),size(X,2));
%   m = chanvese_exp(X,m,150,0.1,nodeBelMF,1,prior_phis, kernel_sigma, Epsilon, beta, lambda); 
%   % mu is standard 0.1
%   phi = make_sdfunc(m);  
%   %pause
% end
% 
% % Plot final result
% matti = reshape(nodeBelMF(:,2),nRows,nCols);
% matti(matti>0.5) = 1;
% matti(matti<=0.5)= 0;
% final_bound = findboundary(matti);
% figure;
% imagesc(X+final_bound);  %imagesc(final_bound);


seg = chanvese_exp(X,m,400,0.1,[],0,prior_phis, 0, Epsilon,beta,lambda); 

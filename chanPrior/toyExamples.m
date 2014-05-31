clear all
close all

addpath('../MRF/');
% addpath('../Chan-Vese/');
load('../data/prior/shapePriorShift.mat');
loadDataFruitfly;
%generateData;

%Run MRF
lambda = 1;
beta = 10;
Epsilon = 1;
prior_phis = Phi(:,:,19:20); 
no_ex = size(prior_phis,3);

% Compute kernel sigma
kernel_sigma = sum(get_Hdistance(circshift(prior_phis,[0 0 -1]),prior_phis,Epsilon));

% initialization
m = zeros(size(X,1),size(X,2));
m(10:size(X,1)-10,10:size(X,2)-10) = 1;
% m0 = prior_phis(:,:,30);
% m(m0>0.5) = 1;
% m(m0<=0.5) = 0;
phi = m(:);
MuEst = [0,1];
SigmaEst = [0.1,0.1];

seg = chanvese(X,m,600,0.1,[],0,prior_phis, kernel_sigma, Epsilon); 

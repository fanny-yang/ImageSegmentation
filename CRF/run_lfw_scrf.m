% This script will train and evaluate spatial conditional random fields
% for the task of Hair/Skin/Background labeling on LFW part label dataset
%
% input (see config_scrf.m for default values)
%   rmposfeat   : remove position features from superpixel features
%   verbose     : display progress during testing
%   dim_crf     : N (see reference)
%   lrl2reg     : weight decay for node weights
%   l2reg_node  : weight decay for node weights
%   l2reg_edge  : weight decay for edge weights
%
% output
%   acc_train   : training accuracy
%   acc_valid   : validation accuracy
%   acc_test    : testing accuracy
%
%
% reference:
% Augmenting CRFs with Boltzmann Machine Shape Priors for Image Labeling, CVPR, 2013.
%

function [acc_train, acc_valid, acc_test] = run_lfw_scrf(rmposfeat,verbose,dim_crf,lrl2reg,l2reg_node,l2reg_edge)

%%% --- default parameter values --- %%%
config_scrf;


%%% --- startup --- %%%
startup;
olddim = 250;   % original LFW image size
nlabel = 3;     % number of segmentation labels

load('sds_large.mat','sds');
load('esds_large.mat','esds');
if rmposfeat,
    sds(65:128) = [];
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% --- spatial logistic regression --- %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('model/slr/');
fname_slr = sprintf('slr_l2r%g_rmposfeat%d_N%d',lrl2reg,rmposfeat,dim_crf);
train_lfw_slr;
save(sprintf('%s/%s.mat',fsave_dir,fname_slr),'w_slr');



% update save directory
fsave_dir = sprintf('%s/%s/',fsave_dir,fname_slr);
if ~exist(fsave_dir,'dir'),
    mkdir(fsave_dir);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% --- spatial conditional random field  --- %%%
%%% --- with mean-field inference         --- %%%
%%% --- train edge weights first, and     --- %%%
%%% --- joint train node and edge weights --- %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('model/scrf/');
fname_scrf = sprintf('scrf_edge_only_dim%d_l2n%g_l2e%g_rmposfeat%d',dim_crf,l2reg_node,l2reg_edge,rmposfeat);
fname_scrf_full = sprintf('scrf_dim%d_l2n%g_l2e%g_rmposfeat%d',dim_crf,l2reg_node,l2reg_edge,rmposfeat);
train_lfw_scrf;
save(sprintf('%s/%s.mat',fsave_dir,fname_scrf),'w_scrf');
if ~skip_nodeup,
    w_scrf = w_scrf_full;
    fname_scrf = fname_scrf_full;
    save(sprintf('%s/%s.mat',fsave_dir,fname_scrf),'w_scrf');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% --- evaluation --- %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n=============================\n');
fprintf('Begin testing! (verbose:%d)\n',verbose);
fprintf('=============================\n\n');

acc_train = eval_lfw_scrf(w_scrf, trainnames, trainnums, sds, esds, verbose);
acc_valid = eval_lfw_scrf(w_scrf, validnames, validnums, sds, esds, verbose);
acc_test = eval_lfw_scrf(w_scrf, testnames, testnums, sds, esds, verbose);

fid = fopen(sprintf('%s/scrf.txt',log_dir),'a+');
fprintf(fid,'acc (val) = %g, acc (test) = %g, (%s)\n',acc_valid,acc_test,fname_scrf);
fclose(fid);

return;

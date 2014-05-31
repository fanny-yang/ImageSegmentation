%% Mean field inference

fprintf('Running Mean Field Inference...\n');
[nodeBelMF,edgeBelMF,logZMF] = UGM_Infer_MeanField(nodePot,edgePot,edgeStruct);

figure;
imagesc(reshape(nodeBelMF(:,2),nRows,nCols));
colormap gray
title('Mean Field Estimates of Marginals');
fprintf('(paused)\n');
pause

% fprintf('Running mean field inference and computing max of marginals\n');
% maxOfMarginalsMFdecode = UGM_Decode_MaxOfMarginals(nodePot,edgePot,edgeStruct,@UGM_Infer_MeanField);
% 
% figure;
% imagesc(reshape(maxOfMarginalsMFdecode,nRows,nCols));
% colormap gray
% title('Max of mean field marginals');
% fprintf('(paused)\n');
% pause


%% Loopy Belief Propagation
% 
% fprintf('Running loopy belief propagation for inference...\n');
% [nodeBelLBP,edgeBelLBP,logZLBP] = UGM_Infer_LBP(nodePot,edgePot,edgeStruct);
% 
% figure;
% imagesc(reshape(nodeBelLBP(:,2),nRows,nCols));
% colormap gray
% title('Loopy Belief Propagation Estimates of Marginals');
% fprintf('(paused)\n');
% pause
% 
% fprintf('Running loopy belief propagation and computing max of marginals\n');
% maxOfMarginalsLBPdecode = UGM_Decode_MaxOfMarginals(nodePot,edgePot,edgeStruct,@UGM_Infer_LBP);
% 
% figure;
% imagesc(reshape(maxOfMarginalsLBPdecode,nRows,nCols));
% colormap gray
% title('Max of Loopy Belief Propagation Marginals');
% fprintf('(paused)\n');
% pause
% 
% fprintf('Running loopy belief propagation for decoding...\n');
% decodeLBP = UGM_Decode_LBP(nodePot,edgePot,edgeStruct);
% 
% figure;
% imagesc(reshape(decodeLBP,nRows,nCols));
% colormap gray
% title('Loopy Belief Propagation Decoding');
% fprintf('(paused)\n');
% pause

function [MuEst, SigmaEst] = normalMixtureEM(MuEst,SigmaEst,X,edgeStruct,probA,iter)
  for i = 1:iter
    [nodePot,edgePot] = producePotentials(X,edgeStruct,MuEst,SigmaEst,probA);
    [nodeBelMF,edgeBelMF,logZMF] = UGM_Infer_MeanField(nodePot,edgePot,edgeStruct);
    [MuEst,SigmaEst] = em_max(nodeBelMF,X(:));
  end

addpath('../MRF/');
addpath('../Chan-Vese/');
loadDataFruitfly;
%generateData;

%Run MRF

edgeStruct = createImgStructure(X,nStates);
lambda = 1;
beta = 10;
   
  % initialization
m = zeros(size(X,1),size(X,2));
m(10:size(X,1)-10,10:size(X,2)-10) = 1;
phi = m(:);
MuEst = [0,1];
SigmaEst = [0.1,0.1];
  
for i = 1:10
  expLambdaPhi = exp(lambda*phi);
  probLogistic = 1./(expLambdaPhi + 1);
  probA = [probLogistic, 1 - probLogistic];
  
  figure;
  imagesc(reshape(probA(:,2),nRows,nCols));
  title('prob(y_i|a_i)');
  colorbar
  %pause

  [MuEst,SigmaEst] = normalMixtureEM(MuEst,SigmaEst,X,edgeStruct,probA,5);
  [nodePot,edgePot] = producePotentials(X,edgeStruct,MuEst,SigmaEst,probA);
  [nodeBelMF,edgeBelMF,logZMF] = UGM_Infer_MeanField(nodePot,edgePot,edgeStruct);  
  
  %figure;
  imagesc(reshape(nodeBelMF(:,2),nRows,nCols));
  title('Mean Field Estimates of Marginals');  
  colorbar
  %pause

  EQ_M = lambda*(probLogistic.*nodeBelMF(:,2) - (1-probLogistic).*nodeBelMF(:,1));
  EQ_M = reshape(EQ_M,size(X,1),size(X,2));
  m = chanvese_hybrid(X,m,100,0.1,'chan',beta*EQ_M); 
  phi = m(:);

  %pause
end

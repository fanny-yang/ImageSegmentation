function [label, model, llh] = fitGMM(X,K,rep)

objVal = zeros(1,rep);
initMat = zeros(size(X,1),K,rep);

for i = 1:rep
  init = kmeanPlusPlusInit(X,K);
  initMat(:,:,i) = init;
  [label, model, llh] = emgm(X,init);
  objVal(i) = llh(end);
end
idxMax = min(find(objVal == max(objVal)));
init = initMat(:,:,idxMax);
[label, model, llh] = emgm(X,init);


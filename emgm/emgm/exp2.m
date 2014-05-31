addpath(genpath('../../'))
labelIdx = 120;

[rowIdx,colIdx] = find(imgLabel == labelIdx);
%figure; scatter(colIdx,rowIdx);hold on;
X = [colIdx,rowIdx]';

rep = 10;
[label, model, llh] = fitGMM(X,2,rep);
figure;scatter(X(1,:),X(2,:),'b'); hold on;
idxTemp1 = find(label == 1);
scatter(X(1,idxTemp1),X(2,idxTemp1),'r');
scatter(model.mu(1,1),model.mu(2,1),'*');
scatter(model.mu(1,2),model.mu(2,2),'o');
if (size(model.mu,2) == 3)
  idxTemp2 = find(label == 2);      
  scatter(X(1,idxTemp2),X(2,idxTemp2),'g');
  scatter(model.mu(1,3),model.mu(2,3),'+');	
end

stop

for i = 1:rep
  init = kmeanPlusPlusInit(X,3);
  [label, model, llh] = emgm(X,init);

  figure; scatter(X(1,:),X(2,:),'b'); hold on;    
  idxTemp1 = find(label == 1);
  scatter(X(1,idxTemp1),X(2,idxTemp1),'r');
  scatter(model.mu(1,1),model.mu(2,1),'*');
  scatter(model.mu(1,2),model.mu(2,2),'o');
  if (size(model.mu,2) == 3)
    idxTemp2 = find(label == 2);      
    scatter(X(1,idxTemp2),X(2,idxTemp2),'g');
    scatter(model.mu(1,3),model.mu(2,3),'+');	
  end
  objVal(i) = llh(end);
title(['LLH=',num2str(llh(end))]);
pause
end


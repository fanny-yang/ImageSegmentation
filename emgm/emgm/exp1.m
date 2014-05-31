addpath(genpath('../../'))
labelIdx = 107;

Itemp = zeros(size(imgLabel));
Itemp(find(imgLabel == labelIdx)) = 1;
imgDist = bwdist(1-Itemp);
figure;imagesc(imgDist);colorbar;

[rowIdx,colIdx] = find(imgLabel == labelIdx);
xlim([min(colIdx),max(colIdx)]);
ylim([min(rowIdx),max(rowIdx)]);
stop
figure; scatter(colIdx,rowIdx);hold on;
X = [colIdx,rowIdx]';
X_demeaned = X - repmat(mean(X,2),1,size(X,2));
[U,S,V] = svd(X_demeaned);

%figure; scatter(X_demeaned(1,:),X_demeaned(2,:));hold on;
D = U;
D(:,1) = D(:,1)*10;
D(:,2) = D(:,2)*10;

%plot([0,D(1,1)],[0,D(2,1)]);
%plot([0,D(1,2)],[0,D(2,2)],'r');

mu1 = D(:,2)+mean(X,2);
mu2 = -D(:,2)+mean(X,2);
%init = [mu1';mu2'];
init = [mu1,mu2];


%init = [485,430;485,450];
%init = [655,260;655,290];
%init = [360,350;130,155];
[label, model, llh] = emgm(X,init);

%idxTemp1 = find(label == 1);
%scatter(colIdx(idxTemp1),rowIdx(idxTemp1),'r');



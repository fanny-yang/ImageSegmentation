shapePrior = imread('../data/prior/shapePrior.png');
shapePrior = mean(shapePrior,3);
shapePrior = shapePrior/255;
shapePrior = imresize(shapePrior,[50,122]);
figure;imagesc(shapePrior);colorbar;
[nRowsPrior, nColsPrior] = size(shapePrior);

idxStartCols = 10:2:20;
idxStartRows = 30:2:40;
Phi = zeros(nRows,nCols,length(idxStartCols)*length(idxStartRows));

k = 1
for i = 1:length(idxStartCols)
  for j = 1:length(idxStartRows)
    idxStartCol = idxStartCols(i);
    idxStartRow = idxStartRows(j);
   
    Y = zeros(nRows,nCols);
    Y(idxStartRow + (1:nRowsPrior), idxStartCol + (1:nColsPrior)) = shapePrior; 
    Phi(:,:,k) = Y;
    k = k + 1;
  end
end

save('../data/prior/shapePirorShift.mat','Phi');

stop
for i = 1:size(Phi,3)
  figure;imagesc(Phi(:,:,i));
  pause;
end
addpath('../utilities');
addpath('../emgm/emgm/');
Img = imread('./insituAT09943_0.jpe');

Img = Img(:,80:end,:);
figure;imshow(Img);
pause

[imgLabel,size_CC] = embryo_watershed(Img);
[nRows,nCols] = size(imgLabel);

idxMax = find(size_CC == max(size_CC));
figure; hist(size_CC(setdiff(1:length(size_CC),idxMax)),100);
pause;

num_CC = length(size_CC);
num_emb = zeros(1,num_CC);

label_count = num_CC+1;
imgLabelOut = imgLabel;
producePlot = 1;
for k = 1:num_CC
  if size_CC(k) < 500
    num_emb(k) = 0;
  elseif size_CC(k) >=500 && size_CC(k)<= 700
    num_emb(k) = 1;
  elseif size_CC(k) >=1000 && size_CC(k)<= 1400
    num_emb(k) = 2;
  elseif size_CC(k) >=1500 && size_CC(k)<= 2100
    num_emb(k) = 3;
  else
    num_emb(k) = 0;
  end

  if num_emb(k) == 0
    continue;
  elseif num_emb(k) == 1
    continue;
  else
    [rowIdx,colIdx] = find(imgLabel == k);			
    X = [colIdx,rowIdx]';
    figure;scatter(X(1,:),X(2,:),'b'); pause;close;
    rep = 10;
    [label, model, llh] = fitGMM(X,num_emb(k),rep);
    
    idxTemp1 = find(label == 1);
    
    colIdxTemp = X(1,idxTemp1);
    rowIdxTemp = X(2,idxTemp1);
    
    idxVecTemp1 = rowIdxTemp + (colIdxTemp-1)*nRows;
    imgLabelOut(idxVecTemp1) = label_count;
   
    imgTemp = zeros(size(imgLabelOut));
    imgTemp(idxVecTemp1) = 1;
    bddTemp = findboundary(imgTemp,2);
    bddIdxTemp = find(bddTemp == 1);
    imgLabelOut(bddIdxTemp) = 0;
    
    label_count = label_count + 1;
    if (size(model.mu,2) == 3)
      idxTemp2 = find(label == 2);      
      colIdxTemp = X(1,idxTemp2);
      rowIdxTemp = X(2,idxTemp2);
      idxVecTemp2 = rowIdxTemp + (colIdxTemp-1)*nRows;
      imgLabelOut(idxVecTemp2) = label_count;
      label_count = label_count + 1;    
      imgTemp = zeros(size(imgLabelOut));
      imgTemp(idxVecTemp2) = 1;
      bddTemp = findboundary(imgTemp,2);
      bddIdxTemp = find(bddTemp == 1);
      imgLabelOut(bddIdxTemp) = 0;
    end

    if producePlot
      figure;scatter(X(1,:),X(2,:),'b'); hold on;
      scatter(X(1,idxTemp1),X(2,idxTemp1),'r');
      scatter(model.mu(1,1),model.mu(2,1),'*');
      scatter(model.mu(1,2),model.mu(2,2),'o');
      if (size(model.mu,2) == 3)
	idxTemp2 = find(label == 2);      
	scatter(X(1,idxTemp2),X(2,idxTemp2),'g');
	scatter(model.mu(1,3),model.mu(2,3),'+');	
      end
    pause
    end
  end

end
labels = unique(imgLabelOut(:));

xCoord = zeros(length(labels),1);
yCoord = xCoord;
for i = 1:length(labels)
  X = zeros(size(imgLabelOut));
  idxTemp = find(imgLabelOut == i);
  X(idxTemp) = 1;
  [xTemp,yTemp] = computeGravityCenter(X);
  xCoord(i) = xTemp;
  yCoord(i) = yTemp;
end

figure;imagesc(imgLabelOut);
text(xCoord,yCoord,num2strBatch(1:length(labels)));

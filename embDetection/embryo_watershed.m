function [imgLabel,size_CC] = embryo_watershed(img)

% Threshold input image and create binary representation with ROI in white
level=graythresh(img);
BW=im2bw(img, level);
BW=~BW;

% Cityblock distance transform
% cf. http://stackoverflow.com/questions/7423596/over-segmentation-in-the-marker-controlled-watershed-in-matlab 
imgDist=-bwdist(~BW,'cityblock');
% Set background pixels to inf
imgDist(~BW)=-inf;

% EF: Get rid of remaining surplus local minima
imgDist2 = imhmin(imgDist, 1);

% Do watershed
imgLabel = watershed(imgDist2);    

labels = unique(imgLabel(:));
xCoord = zeros(length(labels),1);
yCoord = xCoord;
size_CC = zeros(1,length(labels));
for i = 1:length(labels)
  X = zeros(size(BW));
  idxTemp = find(imgLabel == i);
  X(idxTemp) = 1;
  size_CC(i) = length(idxTemp);
  [xTemp,yTemp] = computeGravityCenter(X);
  xCoord(i) = xTemp;
  yCoord(i) = yTemp;
end

% Show results
%figure;imshow(imgLabel==0,'InitialMagnification','fit');
figure;imagesc(imgLabel);
text(xCoord,yCoord,num2strBatch(1:length(labels)));

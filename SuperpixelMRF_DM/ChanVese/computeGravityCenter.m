function [x,y] = computeGravityCenter(phi)
mask = zeros(size(phi));
mask(phi>0) = 1;
mask(phi<=0) = 0;

yDist = sum(mask,2);
yWeight = yDist/sum(yDist);
y = weightMeanVar((1:size(phi,1))',yWeight);

xDist = sum(mask,1);
xWeight = xDist/sum(xDist);
x = weightMeanVar(1:size(phi,2),xWeight);

%figure;imagesc(mask);hold on
%scatter(x,y);
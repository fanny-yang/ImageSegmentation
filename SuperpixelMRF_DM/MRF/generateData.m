clear all;
close all;

Mu = [0,1]; 
Sigma = [0.1,0.4];

Img = imread('../data/shostakovich.jpeg');
Img = imresize(Img,0.25);
Img = mean(Img,3);
Y = (255 - Img)/255;
%Y = Img;
indF = find(Y>.5);
indB = find(Y<=.5);
Y(indF) = 1;
Y(indB) = 0;

noiseF = Sigma(2)*randn(length(indF),1);
noiseB = Sigma(1)*randn(length(indB),1);

X = Y;
X(indF) = Y(indF) + noiseF;
X(indB) = Y(indB) + noiseB;

figure;imshow(1-X);

nNodes = length(X(:));
nStates = 2;
[nRows,nCols] = size(X);




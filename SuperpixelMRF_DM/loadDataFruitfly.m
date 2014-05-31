
doEmbryo = 0;
doGut = 0;
doTouching = 1;

if doEmbryo
  X = imread('../data/embryo.png');
  X = imresize(X,0.5);
elseif doGut
  X = imread('../data/gut.png');
  X = imresize(X,0.5);
  X = mean(X,3);
  %X = histeq(X);
elseif doTouching
  imgName = '../data/touchingEmb/touching2';
  X = imread([imgName,'.jpg']);  
  X = imresize(X,0.3);
  mask = imread('../data/touchingEmb/mask2.png');
  mask = imresize(mask(:,:,1),.3);
end


if superpixel ==0
  X = mean(X,3);
  X = (255 - X)/255;
end

%figure;imagesc(X);
nNodes = length(X(:));
nStates = 2;
[nRows,nCols] = size(X);

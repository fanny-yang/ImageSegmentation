function imgOut = imgZoom(imgIn,r)
if r == 1
  imgOut = imgIn;
elseif r > 1
  [nRows,nCols] = size(imgIn);
  imgTemp = imresize(imgIn,r,'nearest');
  [nRows1,nCols1] = size(imgTemp);
  rowIdx = ceil((nRows1 - nRows)/2+1:(nRows1 + nRows)/2);
  colIdx = ceil((nCols1 - nCols)/2+1:(nCols1 + nCols)/2);
  imgOut = imgTemp(rowIdx,colIdx);
elseif r < 1
  [nRows,nCols] = size(imgIn);
  imgTemp = imresize(imgIn,r,'nearest');
  [nRows1,nCols1] = size(imgTemp);
  rowIdx = ceil((nRows - nRows1)/2+1:(nRows1 + nRows)/2);
  colIdx = ceil((nCols - nCols1)/2+1:(nCols1 + nCols)/2);
  imgOut = zeros(nRows,nCols);
  imgOut(rowIdx,colIdx) = imgTemp;
else
  imgOut = imgIn;
end
 
function I = imgShift(phi,nRowsShift,nColsShift)
% may contain bugs

I = zeros(size(phi));
[nRows,nCols] = size(phi);

if nRowsShift >= 0
  startRowIdx_phi = 1;
  endRowIdx_phi = nRows - nRowsShift;  
  startRowIdx_I = 1 + nRowsShift;
  endRowIdx_I = nRows;
  I(1:(startRowIdx_I-1),:) = repmat(phi(startRowIdx_phi,:),[startRowIdx_I-1,1]);
else
  startRowIdx_I = 1;
  endRowIdx_I = nRows + nRowsShift;
  startRowIdx_phi = 1 - nRowsShift;
  endRowIdx_phi = nRows;  
  I((endRowIdx_I+1):end,:) = repmat(phi(endRowIdx_phi,:),[nRows-endRowIdx_I,1]);
end

if nColsShift >= 0
  startColIdx_phi = 1;
  endColIdx_phi = nCols - nColsShift;  
  startColIdx_I = 1 + nColsShift;
  endColIdx_I = nCols;
  I(:,1:(startColIdx_I-1)) = repmat(phi(:,startColIdx_phi),[1,startColIdx_I-1]);
  
else
  startColIdx_I = 1;
  endColIdx_I = nCols + nColsShift;
  startColIdx_phi = 1 - nColsShift;
  endColIdx_phi = nCols;  
  I(:,(endColIdx_I+1:end)) = repmat(phi(:,endColIdx_phi),[1,nCols - endColIdx_I]);
end

I(startRowIdx_I:endRowIdx_I,startColIdx_I:endColIdx_I) = phi(startRowIdx_phi:endRowIdx_phi,startColIdx_phi:endColIdx_phi);

% input Sp superpixel
addpath('../utilities/');

numPixels = zeros(nSp,1);
se = strel('disk',1);
idxArray = cell(1,nSp);
for i = 1:nSp
  Y = zeros(nRows,nCols);
  idxTemp = find(Sp==i);
  Y(idxTemp) = 1;
  numPixels(i) = length(idxTemp);
  dilatedY = imdilate(Y,se);
  idxArray{i} = find(dilatedY == 1);
end


adj = sparse(nSp,nSp);
for i = 1:(nSp-1)
  for j = (i+1):nSp
    idx1 = idxArray{i};
    idx2 = idxArray{j};
    if length(intersect(idx1,idx2)) > numNbThreshold
      adj(i,j) = 1;
      adj(j,i) = 1;
    end
  end
end

edgeStruct = UGM_makeEdgeStruct(adj,nStates);




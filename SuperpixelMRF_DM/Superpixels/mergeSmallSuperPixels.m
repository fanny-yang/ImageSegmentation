smallSPIdx = find(numPixels <= 100);
for i = 1:length(smallSPIdx)
  idx = smallSPIdx(i);
  nbIdx = find(adj(idx,:)==1);
  nbIdx = setdiff(nbIdx,smallSPIdx);
  overlapTemp = zeros(length(nbIdx),1);
  for j = 1:length(nbIdx)
    nbIdxTemp = nbIdx(j);
    overlapTemp(j) = length(intersect(idxArray{i},idxArray{j}));
  end
  bestNbIdx = min(find(overlapTemp == max(overlapTemp)));
  Sp(find(Sp == idx)) = nbIdx(bestNbIdx);
end

uniqueSp = unique(Sp);
nSp = length(uniqueSp);

for i = 1:nSp
  Sp(find(Sp == uniqueSp(i))) = i;
end

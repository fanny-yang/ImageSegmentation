function outImage = labelSuperpixelsOnImage(Sp,x)
outImage = Sp;
nSp = length(unique(Sp))
for i = 1:nSp
  outImage(find(Sp == i)) = x(i);
end
function boundary_im = findboundary(X,radius)
if ~exist('radius')
            radius = 1;
end
se = strel('disk',radius);
dilatedIm = imdilate(X,se);
boundary_im = dilatedIm - X;

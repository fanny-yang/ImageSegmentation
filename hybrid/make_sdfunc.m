function sdfunc = make_sdfunc(mask)

% Gets you a function with >0 inside the region and <0 outside.
% -> a signed distance function
sdfunc = -(bwdist(mask)-bwdist(1-mask)+im2double(mask)-.5);
function sdfunc = make_sdfunc(mask)
sdfunc = -(bwdist(mask)-bwdist(1-mask)+im2double(mask)-.5);
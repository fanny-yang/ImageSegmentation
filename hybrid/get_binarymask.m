function m = get_binarymask(Phi)
 m = zeros(size(Phi));
 m0 = Phi;
 m(m0>0.5) = 1;
 m(m0<=0.5) = 0;

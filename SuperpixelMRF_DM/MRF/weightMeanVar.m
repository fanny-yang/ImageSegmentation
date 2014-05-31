function [mu, sigma2] = weightMeanVar(x,w)
% x: the image vector
% w: the weights
mu = sum(w.*x)/sum(w);
sigma2 = sum(w.*((x - mu).^2))/sum(w);
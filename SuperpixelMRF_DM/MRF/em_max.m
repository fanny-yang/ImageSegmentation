function [Mu, Sigma] = em_max(dist, X)
% dist: probability distribution matrix, size n by nStates, 
% X: the vectorized image features; n by p
 
[n,p] = size(X);  % number of pixels
l = size(dist,2); % number of labels

Mu = zeros(p,l);
Sigma = zeros(p,l);
for k = 1:l
  w = dist(:,k);
  
  for j = 1:p
  
    [mu,sigma2] = weightMeanVar(X(:,j),w);
    Mu(j,k) = mu;
    Sigma(j,k) = sqrt(sigma2) + eps;
  end
end


% Calculate number of pixels per label
%no_pix = dist*ones(n,l);  % size lx1

% Calculate mean
%mu = dist*img./no_pix; % size lx1

% Calculate variance
%var_vec = (repmat(img,[l 1]) - repmat(mu, [1 n])).^2;  % size lxn
%sigma = dist'*var_vec'./no_pix;


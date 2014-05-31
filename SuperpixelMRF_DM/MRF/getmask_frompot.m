function mask = getmask_frompot(nodePot, nRows, nCols)

% Getting marginal potentials 
% Output: Vector with 1 if prob y =1 > prob y = 0

[val ind] = max(nodePot,[],2); % Get the index with larger prob]
% For ind = 1 set = 1, if ind = 2 then set = 0.

mask = reshape(mod(ind,2),nRows,nCols);
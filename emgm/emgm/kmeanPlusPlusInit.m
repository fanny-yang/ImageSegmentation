function x = kmeanPlusPlusInit(X,K)
% X: the data matrix;
% K: number of data points to be drawn;

if ~exist('method')
    method = 'convex';
end

N = size(X,2);

% choose the first element randomly or according to the centered magnitudes...
sampleIdx = randsample(1:N,1);
restIdx = setdiff(1:N,sampleIdx);
x = X(:,sampleIdx);
xRep = repmat(x,1,N-1);
Xrest = X(:,restIdx);
distTemp = sum((Xrest - xRep).^2);
prob = distTemp/sum(distTemp);

% Now choose the second element...
idx2nd = randsample(restIdx,1,true,prob);
sampleIdx = [sampleIdx,idx2nd];
restIdx = setdiff(1:N,sampleIdx);

% use least squares on prob simplex or K-mean++ to compute sampling
% probabilities for the subsequential archetypes. 
for k = 3:K
    x = X(:,sampleIdx);
    Xrest = X(:,restIdx);
    
    %disMat = slmetric_pw(x,Xrest,'sqdist');
    
    M = bsxfun(@plus, sum(x .* x, 1)', (-2) * x' * Xrest);        
    M = bsxfun(@plus, sum(Xrest .* Xrest, 1), M);        
    M(M < 0) = 0;                        
    disMat = sqrt(M);
    distTemp = min(disMat);
    
    if (sum(distTemp)<0.00001)
        print('All data points have zero distances to the selected...');
        break;
    end
    prob = distTemp/sum(distTemp);
   
    idxknd = randsample(restIdx,1,true,prob);
    sampleIdx = [sampleIdx,idxknd];
    restIdx = setdiff(restIdx,idxknd);
end

x = X(:,sampleIdx);
        


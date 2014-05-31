function [nodePot,edgePot] =  producePotentials(X,edgeStruct,Mu,Sigma,probA)

% X: the vectorized image features; nPixels by p
% Mu: estimate for the means; p by nStates; 
% Sigma: estimate for the standard deviations; p by nStates;
% probA: the estimated prob(y_i|a_i) vector
% create note and edge potentials
% probA = ones(1,nNodes)*(1/nStates);

p = size(X,2);
Xstd = UGM_standardizeCols(X,1);
nNodes = size(X,1);
nStates = size(Mu,2);

nodePot = zeros(nNodes,nStates);

for i = 1:nNodes
  for l = 1:nStates
    expTerm = exp(- 0.5*sum((Sigma(:,l).^(-2)).*(X(i,:)'-Mu(:,l)).^2) + 50);
    nodePot(i,l) = (1/((6.28)^(p/2)*prod(Sigma(:,l))))*expTerm*probA(i,l) + eps;
   %{
    if length(intersect([1,14,33,41],i)) > 0
      [i,l]
      expTerm
      temp = X(i,:)'-Mu(:,l);
      temp
    end
    %}
    %nodePot(i,l) = (1/(sqrt(6.28)*Sigma(l)))*exp(-(Xvec(i)-Mu(l))^2/(2*Sigma(l)^2))*probA(i,l)+eps;
  end
end

edgePot = zeros(nStates,nStates,edgeStruct.nEdges);
for e = 1:edgeStruct.nEdges
  n1 = edgeStruct.edgeEnds(e,1);
  n2 = edgeStruct.edgeEnds(e,2);

  pot_same = exp(1.8 + .3*1/(1+norm(Xstd(n1,:)-Xstd(n2,:))));
  edgePot(:,:,e) = [pot_same 1;1 pot_same];
  %edgePot(:,:,e) = [1 1;1 1];
  
end
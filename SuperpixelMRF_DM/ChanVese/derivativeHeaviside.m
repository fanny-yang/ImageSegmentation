function delta_eps = derivativeHeaviside(phi0,Epsilon,method)

if ~exist('method')
  method = 'sin';
end

if strmatch(method,'sin','exact')
  idxOut = find(abs(phi0)>Epsilon);        
  idxIn = find(abs(phi0)<=Epsilon);        
  delta_eps = zeros(size(phi0));
  delta_eps(idxIn) = (1/(2*Epsilon))*(1+cos(phi0(idxIn)/Epsilon));
else
  delta_eps = (1/Epsilon)*(1./(1+phi0.^2/Epsilon^2));
end         
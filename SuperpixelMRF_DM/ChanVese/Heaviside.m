function H=Heaviside(z,Epsilon,method)

if ~exist('method')
  method = 'sin';
end

%Epsilon=10^(-5);
if strmatch(method,'sin','exact')
  % H=zeros(size(z,1),size(z,2));
  H = zeros(size(z));
  idx1=find(z>Epsilon);
  idx2=find(z<Epsilon & z>-Epsilon);

  H(idx1)=1;
  for i=1:length(idx2)
    H(idx2(i))=1/2*(1+z(idx2(i))/Epsilon+1/pi*sin(pi*z(idx2(i))/Epsilon));
  end;
  
%   idxOut = find(abs(phi0)>Epsilon);        
%   idxIn = find(abs(phi0)<=Epsilon);        
%   delta_eps = zeros(size(phi0));
%   delta_eps(idxIn) = (1/(2*Epsilon))*(1+cos(phi0(idxIn)/Epsilon));
else
  H = .5*(1 + 0.6366*atan(z/Epsilon));
end
function [edgePot,edgePotVec] = produceEdgePotentials(probBoundary,edgeStruct,idxArray,nStates)


edgePot = zeros(nStates,nStates,edgeStruct.nEdges);
edgePotVec = zeros(edgeStruct.nEdges,1);
avgProbBoundary = edgePotVec;
for e = 1:edgeStruct.nEdges
  n1 = edgeStruct.edgeEnds(e,1);
  n2 = edgeStruct.edgeEnds(e,2);
  
       
  idx1 = idxArray{n1};
  idx2 = idxArray{n2};
  idxTemp = intersect(idx1,idx2);
  avgProbBoundary(e) = mean(probBoundary(idxTemp));
%   if (n1 == 47 && n2 == 51)||(n1 == 51 && n2 == 47)
%     avgProbBoundary(e) = 1;
%   end
%   if (n1 == 46 && n2 == 51)||(n1 == 51 && n2 == 46)
%     avgProbBoundary(e) = 1;
%   end
%   if (n1 == 33 && n2 == 5)||(n1 == 5 && n2 == 33)
%     avgProbBoundary(e) = 1;
%   end
%   
end

avgProbBoundary = avgProbBoundary/max(avgProbBoundary);

for e = 1:edgeStruct.nEdges  
  pot_same = 1 - avgProbBoundary(e) + 0.001;
  %pot_same = exp(-15*(avgProbBoundary(e)-0.01));
  %pot_same = 
  edgePotVec(e) = pot_same;
  edgePot(:,:,e) = [pot_same 1;1 pot_same];
  %edgePot(:,:,e) = [1 1;1 1];
end
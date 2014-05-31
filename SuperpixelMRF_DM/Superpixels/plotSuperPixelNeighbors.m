% plot the superpixel index on the segmented image

probThreshold = quantile(edgePotVec,[.3,.6,.9])

spCoord = zeros(nSp,2);
figure; imagesc(I_sp);
for i = 1:nSp
  [y,x] = find(Sp == i);
  spCoord(i,:) = [mean(x),mean(y)];
end

text(spCoord(:,1), spCoord(:,2), num2strBatch(1:nSp)); 
hold on;
for e = 1:edgeStruct.nEdges
  n1 = edgeStruct.edgeEnds(e,1);
  n2 = edgeStruct.edgeEnds(e,2);
  
  if edgePotVec(e)<=probThreshold(1) 
    plot(spCoord([n1,n2],1), spCoord([n1,n2],2),'b');    
  elseif edgePotVec(e)>probThreshold(1) && edgePotVec(e)<=probThreshold(2)
    plot(spCoord([n1,n2],1), spCoord([n1,n2],2),'g');
  elseif edgePotVec(e)>probThreshold(2) && edgePotVec(e)<=probThreshold(3)
    plot(spCoord([n1,n2],1), spCoord([n1,n2],2),'y');
  elseif edgePotVec(e)>probThreshold(3) 
    plot(spCoord([n1,n2],1), spCoord([n1,n2],2),'r');
  end   
 
end    


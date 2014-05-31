function indicator = checkstop_contour(old,new,dt)
% indicate whether we should performance further iteraions or stop

% He set delta  = 0.18
delta = .08;


layer = size(new,3);

for i = 1:layer
    old_{i} = old(:,:,i);
    new_{i} = new(:,:,i);
end

% Get H distance between the two
diff = get_Hdistance(old,new,1e-7); % Basically a real heaviside


% Get difference for the black portion
ind = find(abs(new)<=.5);
M = length(ind);
Q = sum(abs(new(ind)-old(ind)))./M;

% Get Hamming distance
old(old>0)=1;
old(old<=0)=0;
new(new>0)=1;
new(new<=0)=0;
diff_mat = abs(old-new);
ham_dist = sum(sum(diff_mat));


if ham_dist <= 1 && Q<=dt*delta^2 %dt*delta^2
    indicator = 1;
else
    indicator = 0;
end

% if layer
%     ind = find(abs(new)<=.5);
%     M = length(ind);
%     Q = sum(abs(new(ind)-old(ind)))./M;
%     if Q<=dt*delta^2
%         indicator = 1;
%     else
%         indicator = 0;
%     end
% else
%     ind1 = find(abs(old_{1})<1);
%     ind2 = find(abs(old_{2})<1);
%     M1 = length(ind1);
%     M2 = length(ind2);
%     Q1 = sum(abs(new_{1}(ind1)-old_{1}(ind1)))./M1;
%     Q2 = sum(abs(new_{2}(ind2)-old_{2}(ind2)))./M2;
%     if Q1<=dt*delta^2 && Q2<=dt*delta^2
%         indicator = 1;
%     else
%         indicator = 0;
%     end
%end
return

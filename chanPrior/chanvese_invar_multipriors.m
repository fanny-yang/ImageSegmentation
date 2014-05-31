%function seg = chanvese_shiftInvar(I,mask,num_iter,mu, nodeBelMF, hybrid_b, prior_phis, kernel_sigma, Epsilon,beta,lambda)
Epsilon = 3;
kernel_sigma = 1000;
load('../data/prior/embPrior.mat');

%% Prior Pre-processing
Phi = 1 - Phi;
Phi(Phi>0.5) = 1;
Phi(Phi<=0.5) =0;
prior_phis = Phi(:,:,1:3);
num_priors = size(prior_phis,3);

rng(215);
imgTemp = imgZoom(Phi(:,:,1),.9);
imgTemp = imrotate(imgTemp,20,'bilinear','crop');
imgTemp = imgShift(imgTemp,10,10);
%imgTemp(45:50,:) = .4;
imgTemp(:,80:85) = .4;

imgTemp1 = imgZoom(Phi(:,:,1),.9);
imgTemp1 = imrotate(imgTemp1,20,'bilinear','crop');
imgTemp2 = imgZoom(Phi(:,:,2),1);
imgTemp1 = imgShift(imgTemp1,17,0);

imgTemp2 = 1*imgShift(imgTemp2,-21,0);
imgTemp = imgTemp1 + imgTemp2;
imgTemp = imgTemp + 0.2*randn(size(imgTemp));


I = imgTemp;

figure;imagesc(I)

no_ex = size(prior_phis,3);
[m,n] = size(I);


mask = zeros(size(I));
mask(10:size(I,1)-10,5:size(I,2)-10) = 1;
			%	mask = imgShift(imgZoom(prior_phis(:,:,1),.8),20,0);

mask(mask>0) = 1;
mask(mask<=0) = 0;
phi = make_sdfunc(mask);

%   initial force, set to eps to avoid division by zeros
force = eps;

figure(1)

whichHeaviside = 'tan';

A = n/2*ones(1,num_priors); % initial x coord for the priors
B = m/2*ones(1,num_priors); % initial y coord for the prior
R = ones(1,num_priors); % scales...
Theta = zeros(1,num_priors); % initial angle

Psi_0s = zeros(2*m,2*n,num_priors);
Psi_0xs = zeros(2*m,2*n,num_priors);
Psi_0ys = zeros(2*m,2*n,num_priors);

idxRowsTemp = ceil((.5*m+1):(1.5*m));
idxColsTemp = ceil((.5*n+1):(1.5*n));

for i = 1:num_priors
  psi_0 = prior_phis(:,:,i);
  Psi_0 = zeros(2*m,2*n);
  Psi_0(idxRowsTemp,idxColsTemp) = psi_0;
  Psi_0 = make_sdfunc(Psi_0);
  Psi_0s(:,:,i) = Psi_0;
  [Psi_0x, Psi_0y] = gradient(Psi_0);
  Psi_0xs(:,:,i) = Psi_0x;
  Psi_0ys(:,:,i) = Psi_0y;
end


%-- Main loop
num_iter = 2000;
for i=1:num_iter
%{
  seg = phi>0;
  seg = make_sdfunc(phi);
  inidx = intersect(find(seg>0), find(abs(seg)<=Epsilon));
  outidx = intersect(find(seg<=0), find(abs(seg)<=Epsilon));
  activeIdx = union(inidx,outidx);
  c1 = mean(I(inidx));
  c2 = mean(I(outidx));
  force_image = zeros(size(I));
  force_image(activeIdx) = -(I(activeIdx) - c1).^2 + (I(activeIdx) - c2).^2;
  %pause
%}

    inidx = find(phi>0); % frontground index
    outidx = find(phi<=0); % background index
    %force_image = 0; % initial image force for each layer
    
    %% Image dependent force p(i|phi)
   
    L = im2double(I);
    c1 = mean(I(inidx));
    c2 = mean(I(outidx));
    %c1 = sum(sum(L.*Heaviside(phi,Epsilon,whichHeaviside)))/(length(inidx)+eps); % average inside of Phi0
    %c2 = sum(sum(L.*(1-Heaviside(phi,Epsilon,whichHeaviside))))/(length(outidx)+eps); % average outside of Phi0
    force_image=-(L-c1).^2+(L-c2).^2;
    %figure(2);imagesc(force_image);colorbar;
    %pause

    delta_A = A*0.0;
    delta_B = B*0.0;
    delta_Theta = Theta*0.0;
    delta_R = R*0.0;
   
    xMat = repmat(1:size(phi,2),[size(phi,1),1]);
    yMat = repmat((1:size(phi,1))',[1,size(phi,2)]);
    
    prior_term = 0;
    delta_Psi = zeros(m,n,num_priors);
    Psi = zeros(m,n,num_priors);
    Alpha = zeros(1,num_priors);
    for j = 1:num_priors
      a = A(j);
      b = B(j);
      theta = Theta(j);
      r = R(j);

      Psi_0 = Psi_0s(:,:,j);
      Psi_0x = Psi_0xs(:,:,j);
      Psi_0y = Psi_0ys(:,:,j);
      
      xShift = (ceil(a - n/2 + eps) - 1);
      yShift = (ceil(b - m/2 + eps) - 1);
    
      psi = imgZoom(Psi_0,r);
      psi = imrotate(psi,-57.296*theta,'bilinear','crop');
      psi = imgShift(psi,yShift, xShift);
      psi = psi(idxRowsTemp,idxColsTemp);
      Psi(:,:,j) = psi;

      psi_0x = imgZoom(Psi_0x,r);
      psi_0x = imrotate(psi_0x,-57.296*theta,'bilinear','crop');
      psi_0x = imgShift(psi_0x,yShift, xShift);
      psi_0x = psi_0x(idxRowsTemp,idxColsTemp);

      psi_0y = imgZoom(Psi_0y,r);
      psi_0y = imrotate(psi_0y,-57.296*theta,'bilinear','crop');
      psi_0y = imgShift(psi_0y,yShift, xShift);
      psi_0y = psi_0y(idxRowsTemp,idxColsTemp);
          
      %% Compute Shape prior term

      H_phi = Heaviside(phi, Epsilon, whichHeaviside);
      H_psi = Heaviside(psi, Epsilon, whichHeaviside);
      diff_phi_mat = H_phi - H_psi;
      delta_psi = derivativeHeaviside(psi,Epsilon,whichHeaviside);
      delta_phi = derivativeHeaviside(phi,Epsilon,whichHeaviside); %% outside!
      delta_Psi(:,:,j) = delta_psi;

    % compute da/dt db/dt
      delta_a_temp = -(diff_phi_mat.*(cos(theta)*psi_0x - sin(theta)*psi_0y).*delta_psi);
      delta_b_temp = -(diff_phi_mat.*(sin(theta)*psi_0x + cos(theta)*psi_0y).*delta_psi);
    
      delta_a = sum(delta_a_temp(:));
      delta_b = sum(delta_b_temp(:));
     
				% compute d theta/ dt         
      
      xStar = (xMat - xShift)*cos(theta) + (yMat - yShift)*sin(theta);
      yStar = -(xMat - xShift)*sin(theta) + (yMat - yShift)*cos(theta);
    

      delta_theta_temp = -(diff_phi_mat.*(-psi_0x.*yStar + psi_0y.*xStar).*delta_psi);
      delta_theta = sum(delta_theta_temp(:));

    % compute dr /dt
      delta_r_temp = -(diff_phi_mat.*(- psi + psi_0x.*xStar + psi_0y.*yStar).*delta_psi);
      delta_r = sum(delta_r_temp(:));

      alpha_i = exp((- sum(diff_phi_mat(:).^2))/(2*kernel_sigma));
      Alpha(j) = alpha_i;
      %alpha_i = 1;
      prior_term = prior_term + alpha_i*(diff_phi_mat.*delta_phi);

      delta_A(j) = delta_a;
      delta_B(j) = delta_b;
      delta_Theta(j) = delta_theta;
      delta_R(j) = delta_r;
 
    end
        % Normalize prior term if big enough
        if max(max(prior_term)) > 1e-10
            prior_term = prior_term/max(max(prior_term));
        end
        
        if i==1
            
            %-- Display settings
            figure(1);
            subplot(2,2,1); imshow(I); title('Input Image');
            subplot(2,2,2); contour(flipud(phi), [0 0], 'r','LineWidth',1); title('initial contour');
            subplot(2,2,3); imshow(I); title('Segmentation');
            %-- End Display original image and mask
            
        end
    
    
  
    %% Add to external force p(phi)
    
    % Stepsize
    lambda = 1;
    mu = 1;
    
    %delta_phi2 = derivativeHeaviside(phi,10,whichHeaviside); %% outside!

    temp = kappa(phi, Epsilon);
    force1 = mu*temp./max(max(abs(temp))) + force_image;
	%force1 = 0;
    force =  force1 - lambda*prior_term; 
	%figure; imagesc(force1); colorbar; title('force1');
	%figure; imagesc(prior_term); colorbar; title('prior');
	%pause
    
    
    % Normalize the force
    force = force./(max(max(abs(force))));
    
    % Stepsize dt
    dt=0.5;
    
    % Update phi
    old = phi;
    phi = phi+dt.*force;
    new = phi;

    dt1 = 0.01; % magic step size
    A = A +dt1.*delta_A;
    B = B +dt1.*delta_B;
    %pause
    
    dt2 = 0.000005;    
    Theta = mod(double(Theta +dt2*delta_Theta),6.28);

    dt3 = 0.00001;
    R = R + dt3*delta_R;
    R(find(R>1.1)) = 1.1;
    R(find(R<.9)) = .9;

				%theta = pi/3;
    %pause
    %theta =0;

    
    
    
    % Check stopping condition
    ind = find(abs(old)<=Epsilon);
    diff_old_new = sum(abs(old(ind) - new(ind)))/length(ind);
    if (diff_old_new <= 0.002*dt)
      indicator = 1;
    else
      indicator = 0;
    end

    current_cont = zeros(size(phi));
    current_cont(phi<0) = 0;
    current_cont(phi>=0) = 1;
    current_bound = findboundary(current_cont);

				%pause
    
    % Intermediate output
    if (mod(i,20) == 0)
      figure(1)
      for j = 1:num_priors
        subplot(2,2,j);
	imagesc(delta_Psi(:,:,j)); 
	mainTitle = ['theta=',num2str(mod(180*Theta(j)/pi,360)),' (degree); scale=', num2str(R(j))];
	title(mainTitle);
      end
      subplot(2,2,4);
      %showphi(I,phi,i);
      bestPriorIdx = find(Alpha == max(Alpha));
      bestPrior = Psi(:,:,bestPriorIdx);
      idxTemp = intersect(find(bestPrior>=-Epsilon),find(phi>=0)); 
      %idxTemp = find(phi>=0); 
      seg = zeros(m,n);
      seg(idxTemp) = 1;
      bdd = findboundary(seg);
      imagesc(I+2*bdd);
      mainTitle = ['iter=',num2str(i),'; diff = ', num2str(diff_old_new)];
      title(mainTitle);
    end;
    
    if indicator % decide to stop or continue
      break      
    end
   
end

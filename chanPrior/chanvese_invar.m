function [seg,phi,param] = chanvese_invar(I,mask,num_iter,mu, nodeBelMF_CV, prior_phis, kernel_sigma, Epsilon, beta,lambda,priorWeight,param)

num_priors = size(prior_phis,3);
[m,n] = size(I);

mask(mask>0) = 1;
mask(mask<=0) = 0;
phi = make_sdfunc(mask);

%   initial force, set to eps to avoid division by zeros
%force = eps;

figure(1)

whichHeaviside = 'tan';

if ~exist('param')
  A = n/2*ones(1,num_priors); % initial x coord for the priors
  B = m/2*ones(1,num_priors); % initial y coord for the prior
  R = ones(1,num_priors); % scales...
  Theta = zeros(1,num_priors); % initial angle
else
  A = param.A;
  B = param.B;
  R = param.R;
  Theta = param.Theta;
end

Psi_0s = zeros(4*m,4*n,num_priors);
Psi_0xs = zeros(4*m,4*n,num_priors);
Psi_0ys = zeros(4*m,4*n,num_priors);

idxRowsTemp = ceil((1.5*m+1):(2.5*m));
idxColsTemp = ceil((1.5*n+1):(2.5*n));

for i = 1:num_priors
  psi_0 = prior_phis(:,:,i);
  Psi_0 = zeros(4*m,4*n);
  Psi_0(idxRowsTemp,idxColsTemp) = psi_0;
  Psi_0 = make_sdfunc(Psi_0);
  Psi_0s(:,:,i) = Psi_0;
  [Psi_0x, Psi_0y] = gradient(Psi_0);
  Psi_0xs(:,:,i) = Psi_0x;
  Psi_0ys(:,:,i) = Psi_0y;
end


%-- Main loop

for i=1:num_iter
    
    inidx = find(phi>0); % frontground index
    outidx = find(phi<=0); % background index
    
    c1 = mean(I(inidx));
    c2 = mean(I(outidx));
    force_image=-(I-c1).^2+(I-c2).^2;

        % Calculate EQ_M
    if beta > 0 
      expLambdaPhi = exp(lambda*phi(:));
      probLogistic = 1./(expLambdaPhi + 1);
      EQ_M = lambda*(probLogistic.*nodeBelMF_CV(:,2) - (1-probLogistic).*nodeBelMF_CV(:,1));
      EQ_M = reshape(EQ_M,m,n);
    else
      EQ_M = 0;
    end

    %force_image = 0; % initial image force for each layer    
    %L = im2double(I);
    %c1 = sum(sum(L.*Heaviside(phi,Epsilon,whichHeaviside)))/(length(inidx)+eps); % average inside of Phi0
    %c2 = sum(sum(L.*(1-Heaviside(phi,Epsilon,whichHeaviside))))/(length(outidx)+eps); % average outside of Phi0
    %force_image=-(L-c1).^2+(L-c2).^2+force_image;


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
      delta_phi = derivativeHeaviside(phi,Epsilon,whichHeaviside);
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
            subplot(3,2,1); imshow(I); title('Input Image');
            subplot(3,2,2); contour(flipud(phi), [0 0], 'r','LineWidth',1); title('initial contour');
            subplot(3,2,3); imshow(I); title('Segmentation');
            %-- End Display original image and mask
            
        end
    
    
  
    %% Add to external force p(phi)
    
    % Stepsize
    %priorWeight = 0.5;
    
    temp = kappa(phi, Epsilon);
    force1 = mu*temp./max(max(abs(temp))) + force_image;
    force =  force1 + beta*EQ_M - priorWeight*prior_term; 
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
    %dt1 = 0;
    A = A +dt1.*delta_A;
    B = B +dt1.*delta_B;
    
    dt2 = 0.000005;
    %dt2 = 0;
    Theta = mod(double(Theta +dt2*delta_Theta),6.28);

    dt3 = 0.00001;
    %dt3 = 0; %%!!
    R = R + dt3*delta_R;
    R(find(R>1)) = 1;
    R(find(R<.9)) = 0.9;
    
    % Check stopping condition
    ind = find(abs(old)<=Epsilon);
    diff_old_new = sum(abs(old(ind) - new(ind)))/(length(ind)+eps);
    if (diff_old_new <= 0.002*dt)
      indicator = 1;
    else
      indicator = 0;
    end

    indicator = 0;
    % Intermediate output
    if (mod(i,20) == 0)
      figure(1)
      for j = 1:num_priors
        subplot(3,2,j);
	imagesc(delta_Psi(:,:,j)); 
	mainTitle = ['theta=',num2str(mod(180*Theta(j)/pi,360)),' (degree); scale=', num2str(R(j))];
	title(mainTitle);
      end
      
      subplot(3,2,4);
      bestPriorIdx = find(Alpha == max(Alpha));
      bestPrior = Psi(:,:,bestPriorIdx);
      %idxTemp = find(phi>=0); 
      idxTemp = intersect(find(bestPrior>=-Epsilon),find(phi>=0)); 
      seg = zeros(m,n);
      seg(idxTemp) = 1;
      bdd = findboundary(seg);
      imagesc(I+2*bdd);
      mainTitle = ['iter=',num2str(i),'; diff = ', num2str(diff_old_new)];
      title(mainTitle);
      subplot(3,2,5);
      idxTemp = find(phi>=0); 
      seg = zeros(m,n);
      seg(idxTemp) = 1;
      bdd = findboundary(seg);
      imagesc(I+2*bdd);
            
    end;
    
    if indicator 
      break      
    end
   
end

param.A = A;
param.B = B;
param.Theta = Theta;
param.R = R;


bestPriorIdx = find(Alpha == max(Alpha));
bestPrior = Psi(:,:,bestPriorIdx);
idxTemp = intersect(find(bestPrior>=-Epsilon),find(phi>=0)); 
seg = zeros(m,n);
seg(idxTemp) = 1;

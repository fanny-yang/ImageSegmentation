%function seg = chanvese_shiftInvar(I,mask,num_iter,mu, nodeBelMF, hybrid_b, prior_phis, kernel_sigma, Epsilon,beta,lambda)
Epsilon = 3;
load('../data/prior/embPrior.mat');

%% Prior Pre-processing
Phi = 1 - Phi;
Phi(Phi>0.5) = 1;
Phi(Phi<=0.5) =0;
prior_phis = Phi(:,:,:);

rng(215);
imgTemp = imgZoom(prior_phis(:,:,1),.9);
imgTemp = imrotate(imgTemp,20,'bilinear','crop');
imgTemp = imgShift(imgTemp,10,10);
%imgTemp(45:50,:) = .4;
imgTemp(:,80:85) = .4;

imgTemp = imgTemp + .1*randn(size(imgTemp));
I = imgTemp;
figure;imagesc(I)

no_ex = size(prior_phis,3);
[m,n] = size(I);


mask = zeros(size(I));
mask(10:size(I,1)-10,5:size(I,2)-10) = 1;

mask(mask>0) = 1;
mask(mask<=0) = 0;
phi = make_sdfunc(mask);

%   initial force, set to eps to avoid division by zeros
force = eps;

figure(1)

whichHeaviside = 'tan';

a = n/2; % initial x coord for the prior
b = m/2; % initial y coord for the prior
r = 1; % scales...
theta = double(0); % initial angle

psi_0 = prior_phis(:,:,1);
Psi_0 = zeros(2*size(psi_0));


idxRowsTemp = ceil((.5*m+1):(1.5*m));
idxColsTemp = ceil((.5*n+1):(1.5*n));

Psi_0(idxRowsTemp,idxColsTemp) = psi_0;
Psi_0 = make_sdfunc(Psi_0);
[Psi_0x, Psi_0y] = gradient(Psi_0);



%-- Main loop
num_iter = 2000;
for i=1:num_iter
    
    inidx = find(phi>0); % frontground index
    outidx = find(phi<=0); % background index
    force_image = 0; % initial image force for each layer
    
    %% Image dependent force p(i|phi)
   
    L = im2double(I);
    c1 = sum(sum(L.*Heaviside(phi,Epsilon,whichHeaviside)))/(length(inidx)+eps); % average inside of Phi0
    c2 = sum(sum(L.*(1-Heaviside(phi,Epsilon,whichHeaviside))))/(length(outidx)+eps); % average outside of Phi0
    force_image=-(L-c1).^2+(L-c2).^2+force_image;
    

    xShift = (ceil(a - n/2 + eps) - 1)
    yShift = (ceil(b - m/2 + eps) - 1)
    

    psi = imgZoom(Psi_0,r);
    psi = imrotate(psi,-57.296*theta,'bilinear','crop');
    psi = imgShift(psi,yShift, xShift);
    psi = psi(idxRowsTemp,idxColsTemp);
    psiTemp = psi;

    psi_0x = imgZoom(Psi_0x,r);
    psi_0x = imrotate(psi_0x,-57.296*theta,'bilinear','crop');
    psi_0x = imgShift(psi_0x,yShift, xShift);
    psi_0x = psi_0x(idxRowsTemp,idxColsTemp);
    

    psi_0y = imgZoom(Psi_0y,r);
    psi_0y = imrotate(psi_0y,-57.296*theta,'bilinear','crop');
    psi_0y = imgShift(psi_0y,yShift, xShift);
    psi_0y = psi_0y(idxRowsTemp,idxColsTemp);
    
    
      
    %% Compute Shape prior term
    % only one prior psi (psi_0)
    H_phi = Heaviside(phi, Epsilon, whichHeaviside);
    H_psi = Heaviside(psi, Epsilon, whichHeaviside);
    diff_phi_mat = H_phi - H_psi;
    delta_psi = derivativeHeaviside(psi,Epsilon,whichHeaviside);
    delta_phi = derivativeHeaviside(phi,Epsilon,whichHeaviside);
    
    % compute da/dt db/dt
    delta_a_temp = -(diff_phi_mat.*(cos(theta)*psi_0x - sin(theta)*psi_0y).*delta_psi);
    delta_b_temp = -(diff_phi_mat.*(sin(theta)*psi_0x + cos(theta)*psi_0y).*delta_psi);
    
    delta_a = sum(delta_a_temp(:));
    delta_b = sum(delta_b_temp(:));
    

				% compute d theta/ dt    
    xMat = repmat(1:size(phi,2),[size(phi,1),1]);
    yMat = repmat((1:size(phi,1))',[1,size(phi,2)]);
    
    %xStar = (xMat - a)*cos(theta) + (yMat - b)*sin(theta);
    %yStar = -(xMat - a)*sin(theta) + (yMat - b)*cos(theta);
    
    xStar = (xMat - xShift)*cos(theta) + (yMat - yShift)*sin(theta);
    yStar = -(xMat - xShift)*sin(theta) + (yMat - yShift)*cos(theta);
    

    delta_theta_temp = -(diff_phi_mat.*(-psi_0x.*yStar + psi_0y.*xStar).*delta_psi);
    delta_theta = sum(delta_theta_temp(:));

    % compute dr /dt
    delta_r_temp = -(diff_phi_mat.*(- psi + psi_0x.*xStar + psi_0y.*yStar).*delta_psi);
    delta_r = sum(delta_r_temp(:));

 %{
    figure; imagesc(diff_phi_mat);colorbar;
    figure; imagesc(delta_theta_temp);colorbar;
    figure; imagesc(psi);colorbar;
    pause
%}
    prior_term = diff_phi_mat.*delta_phi;
 

   %{
	figure;imagesc(shiftInvarTerm);colorbar
        %figure;imagesc(sum(alpha_i.*diff_phi_mat,3));colorbar
        figure;imagesc(sum(alpha_i.*diff_phi_mat,3)+shiftInvarTerm);colorbar
	temp = phi0;
	temp(temp>0) = 1;
	temp(temp<=0) = 0;
	bddTemp = findboundary(temp);
	figure; imagesc(phi0/max(phi0(:))+bddTemp); 
	pause
        %}

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
    mu = .3;

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
    a = a +dt1.*delta_a;
    b = b +dt1.*delta_b;
    %pause
    
    dt2 = 0.00002;    
    theta = mod(double(theta +dt2*delta_theta),6.28);

    dt3 = 0.00001;
    r = r + dt3*delta_r;
    
				%theta = pi/3;
    %pause
    %theta =0;

    
    
    
    % Check stopping condition
    %indicator = checkstop(old,new,dt);
    indicator = 0;
    current_cont = zeros(size(phi));
    current_cont(phi<0) = 0;
    current_cont(phi>=0) = 1;
    current_bound = findboundary(current_cont);

    %{
    figure(2)
    subplot(2,2,1); imagesc(prior_term+current_bound);colorbar;
    subplot(2,2,2); imagesc(force+current_bound); colorbar;
    subplot(2,2,3); imagesc(phi0/max(max(phi0))+current_bound);   colorbar;
    subplot(2,2,4); imagesc(EQ_M/max(max(EQ_M))+current_bound); colorbar;
    %}
				%pause
    
    % Intermediate output
    if(mod(i,20) == 0)
        figure(1)
        subplot(2,2,3);
        showphi(I,phi,i);
        subplot(2,2,1); imagesc(phi/max(max(phi)));%+imboundary);
	subplot(2,2,4); imagesc(delta_psi); 
	mainTitle = ['theta=',num2str(mod(180*theta/pi,360)),' (degree); scale=', num2str(r)];
	title(mainTitle);
        
    end;
    if indicator % decide to stop or continue
      figure(1)
      showphi(I,phi,i);
      
				% Get mask from level set function phi
      seg = phi>=0; % !!
      figure(1)
      subplot(2,2,4); imshow(seg); title('Global Region-Based Segmentation');
      
    end
   
end


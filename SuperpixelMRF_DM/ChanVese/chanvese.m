function seg = chanvese(I,mask,num_iter,mu, EQ_M, hybrid_b, prior_phis, kernel_sigma, Epsilon)

no_ex = size(prior_phis,3);
method = 'chan';
[m,n] = size(I);

%%  Preprocessing
% Resize
s = 200./min(size(I,1),size(I,2)); % resize scale
if s<1
    I = imresize(I,s);
end

% Some more preprocessing
%if size(I,3)== 3
%    P = rgb2gray(uint8(I));
%    P = double(P);
%elseif size(I,3) == 2
%    P = 0.5.*(double(I(:,:,1))+double(I(:,:,2)));
%else
P = double(I);
%end
layer = 1;

%%  Different masks
if ischar(mask)
    switch lower (mask)
        case 'small'
            mask = maskcircle2(I,'small');
        case 'medium'
            mask = maskcircle2(I,'medium');
        case 'large'
            mask = maskcircle2(I,'large');
        case 'whole'
            mask = maskcircle2(I,'whole');
            %mask = init_mask(I,30);
        case 'whole+small'
            m1 = maskcircle2(I,'whole');
            m2 = maskcircle2(I,'small');
            mask = zeros(size(I,1),size(I,2),2);
            mask(:,:,1) = m1(:,:,1);
            mask(:,:,2) = m2(:,:,2);
        otherwise
            error('unrecognized mask shape name (MASK).');
    end
else
    if s<1
        mask = imresize(mask,s);
    end
    if size(mask,1)>size(I,1) || size(mask,2)>size(I,2)
        error('dimensions of mask unmatch those of the image.')
    end
    %switch lower(method)
    %    case 'multiphase'
    %        if  (size(mask,3) == 1)
    %            error('multiphase requires two masks but only gets one.')
    %        end
    %end
    
end

%% Run algorithm

%-- Initialization
%   Get the level set function of the initial mask
mask = mask(:,:,1);
phi0 = -(bwdist(mask)-bwdist(1-mask)+im2double(mask)-.5);
%   initial force, set to eps to avoid division by zeros
force = eps;
%-- End Initialization

%-- Display settings
figure();
subplot(2,2,1); imshow(I); title('Input Image');
subplot(2,2,2); contour(flipud(phi0), [0 0], 'r','LineWidth',1); title('initial contour');
subplot(2,2,3); title('Segmentation');
%-- End Display original image and mask

%-- Main loop
for i=1:num_iter
    
    inidx = find(phi0>=0); % frontground index
    outidx = find(phi0<0); % background index
    force_image = 0; % initial image force for each layer
    
    %% Image dependent force p(i|phi)
    % Only relevant if i is actually the intensity
    
    for j=1:layer
        L = im2double(P(:,:,j)); % get one image component
        c1 = sum(sum(L.*Heaviside(phi0,Epsilon)))/(length(inidx)+eps); % average inside of Phi0
        c2 = sum(sum(L.*(1-Heaviside(phi0,Epsilon))))/(length(outidx)+eps); % average outside of Phi0
        force_image=-(L-c1).^2+(L-c2).^2+force_image;
        % sum Image Force on all components (used for vector image)
        % if 'chan' is applied, this loop become one sigle code as a
        % result of layer = 1
    end
    
    %% Compute Shape prior term

    
    if kernel_sigma         % Want to incorporate prior?
    
    diff_phi_mat = Heaviside(repmat(phi0, [1 1 no_ex]), Epsilon) - Heaviside(prior_phis, Epsilon);
    %figure; imagesc(diff_phi_mat(:,:,2)); colorbar;
    %figure; imagesc(phi0); colorbar;
    %pause
    % Calculate difference
    delta_eps = (1/Epsilon)*(1./(1+phi0.^2/Epsilon^2));
    no_ex = size(prior_phis,3);

    % Calculate distance
    dist_phis = get_Hdistance(repmat(phi0, [1 1 no_ex]),prior_phis, Epsilon);
    help_vec = exp(-dist_phis/(2*kernel_sigma));
    norm_const = sum(help_vec);
    
    % Actual term to add
    alpha_i = reshape(kron(help_vec,ones(m,n)),[m n no_ex]);
    prior_term = delta_eps.*(sum(alpha_i.*diff_phi_mat,3))/(norm_const*2*kernel_sigma);
    %prior_term = (sum(alpha_i.*diff_phi_mat,3))/(norm_const*2*kernel_sigma);
    
    end
    
    
    %% Add to external force p(phi)
    if hybrid_b ==0 && kernel_sigma == 0         % Want to incorporate MRF?
        force = mu*kappa(phi0, Epsilon)./max(max(abs(kappa(phi0, Epsilon))))+1/layer.*force_image;
    elseif hybrid_b == 0 && kernel_sigma ~= 0
        
        %force1 = mu*kappa(phi0, Epsilon)./max(max(abs(kappa(phi0, Epsilon))))+1e2/layer.*force_image; 
        force1=0;
        force = force1 + 3*(1e5)*prior_term;  
        %figure;imagesc(force1);colorbar;
        %figure;imagesc(force);colorbar;
        %figure;imagesc(prior_term);colorbar;
        %pause
        %pause
        %stopcrap
    else
        force = mu*kappa(phi0, Epsilon)./max(max(abs(kappa(phi0, Epsilon))))+1/layer.*force_image +EQ_M - prior_term;  
    end
    
    
    % Normalize the force
    force = force./(max(max(abs(force))));
    
    % Stepsize dt
    dt=0.5;
    
    % Update phi 
    old = phi0;
    phi0 = phi0+dt.*force;
    new = phi0;
    
    % Check stopping condition
    indicator = checkstop_contour(old,new,dt);
    %indicator = 0;
    
    % Intermediate output
    if(mod(i,20) == 0)
        showphi(I,phi0,i);
        
    end;
    if indicator % decide to stop or continue
        showphi(I,phi0,i);
        
        % Get mask from level set function phi
        seg = phi0>=0; % !!
        
        subplot(2,2,4); imshow(seg); title('Global Region-Based Segmentation');
        
        return;
    end
end;
showphi(I,phi0,i);

% Get mask from level set function phi
seg = phi0<=0; 

subplot(2,2,4); imshow(seg); title('Global Region-Based Segmentation');

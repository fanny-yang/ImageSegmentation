function Hdist = get_Hdistance(phi1_mat,phi2_mat, Epsilon)

if size(phi1_mat,3) ~= size(phi2_mat,3)
    disp('Warning: Dimensions are not equal to each other');
else
    no_phi = size(phi1_mat,3);
    diff_phi_mat = Heaviside(phi1_mat, Epsilon) - Heaviside(phi2_mat, Epsilon);
    
    % Calculate distance for alpha_i
    dist_phis = sum(sum(diff_phi_mat.^2,1),2);
    
    % Reshape to vector
    Hdist = zeros(1,no_phi);
    Hdist(:) = dist_phis(1,1,:);
end


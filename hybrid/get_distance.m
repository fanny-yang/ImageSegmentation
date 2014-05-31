function dist_phis = get_Hdistance(phi1_mat,phi2_mat, Epsilon)


diff_phi_mat = Heaviside(phi1_mat, Epsilon) - Heaviside(phi2_mat, Epsilon);
    
% Calculate distance for alpha_i
dist_phis = sum(sum(diff_phi_mat.^2,1),2); 

function [pi_semiz,N_semiz] = create_pi_semiz(pi_l,pi_g,pi_z3,n_s,n_l,n_g,n_z3,N_j,N_i)

if ~isequal(size(pi_l),[n_l,n_l,n_s,N_j,N_i])
    error('Size pi_l is not (l,l'',s,j,ptype) ')
end
if ~isequal(size(pi_g),[n_g,n_g,n_l])
    error('Size pi_g is not (g,l'',l) ')
end
if ~isequal(size(pi_z3),[n_z3,n_z3])
    error('Size pi_z3 is not (n_z3,n_z3) ')
end

N_semiz = n_l*n_g*n_z3;
pi_semiz = zeros(N_semiz,N_semiz,n_s,N_j,N_i);

for i_c=1:N_i
for j_c=1:N_j
for s_c=1:n_s
    for z_c = 1:N_semiz
        [l_c,g_c,z3_c] = ind2sub([n_l,n_g,n_z3],z_c);
        for zp_c = 1:N_semiz
            [lp_c,gp_c,z3p_c] = ind2sub([n_l,n_g,n_z3],zp_c);
            pi_semiz(z_c,zp_c,s_c,j_c,i_c) = pi_l(l_c,lp_c,s_c,j_c,i_c)*...
                pi_g(g_c,gp_c,l_c)*pi_z3(z3_c,z3p_c);
        end
    end
end
end
end

% Check rows sums to one
for ii=1:N_i
    for j_c=1:N_j
    for s_c=1:n_s
        prob = pi_semiz(:,:,s_c,j_c,ii);
        if any(abs(sum(prob,2)-1)>1e-10)
            warning('Rows of pi_semiz do not sum to 1 for s=%d, age=%d and type=%d \n',s_c,j_c,ii)
        end
    end
    end
end

end %end function

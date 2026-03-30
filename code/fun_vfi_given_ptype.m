function [V,pol_s_idx,pol_s,pol_aprime,pol_aprime_val] = fun_vfi_given_ptype(n_a,n_s,n_l,n_g,N_j,...
    a_grid,s_grid,l_grid,g_grid,pi_l,pi_g,Params)

V          = zeros(n_a,n_l,n_g,N_j+1);
pol_s_idx  = zeros(n_a,n_l,n_g,N_j);
pol_s      = zeros(n_a,n_l,n_g,N_j);
pol_aprime = zeros(n_a,n_l,n_g,N_j);

% Precompute asset grids in column/row form so that f_Return_cpu returns
% a matrix with rows indexed by a' and columns indexed by current assets a.
aprime_col = a_grid;
a_row      = a_grid';

% Backward induction over age.
for j_c = N_j:-1:1
    fprintf('Age = %d \n',j_c)
    V_next = V(:,:,:,j_c+1); % (a',l',g')
    V_s          = -Inf(n_a,n_l,n_g,n_s);
    pol_aprime_s = zeros(n_a,n_l,n_g,n_s);

    % Conditional on each search choice, solve the inner maximization over a'.
    for s_c=1:n_s
        s_val  = s_grid(s_c);
        for l_c=1:n_l
            l_val = l_grid(l_c);
            for g_c=1:n_g
                g_val = g_grid(g_c);
                % Employment transitions depend on search effort and on age.
                pi_l_vec = pi_l(l_c,:,s_c,j_c);
                % Skill transitions depend on current employment.
                pi_g_vec = pi_g(g_c,:,l_c);
                % Expected continuation value as a function of a'.
                EV = f_expval(V_next,pi_l_vec,pi_g_vec); % (a',1)
                % Current-period return matrix indexed by (a',a).
                RetMat_s_z = f_Return_cpu(s_val,aprime_col,a_row,l_val,g_val,...
                    Params.agej(j_c),Params.educ_i,Params.Jr,Params.r,Params.w,Params.ben,Params.pens,Params.gamma,Params.B_s);
                % Bellman RHS for the current (l,g,s) block.
                entireRHS = RetMat_s_z + Params.beta*EV; % (a',a)
                [Vtemp,maxindex]   = max(entireRHS,[],1);  % (1,a)
                % Value and a' policy conditional on the current search choice.
                V_s(:,l_c,g_c,s_c) = Vtemp;
                pol_aprime_s(:,l_c,g_c,s_c) = maxindex;
            end %end g
        end %end l
    end %end s

    % For the current age j, maximize over search effort.
    [V_jj,max_s] = max(V_s,[],4); % (a,l,g)
    V(:,:,:,j_c) = V_jj;
    pol_s_idx(:,:,:,j_c) = max_s;
    pol_s(:,:,:,j_c) = s_grid(max_s);

    % Recover the optimal a' index using the optimal search choice max_s.
    for a_c=1:n_a 
        for l_c=1:n_l
            for g_c=1:n_g
                pol_aprime(a_c,l_c,g_c,j_c) = pol_aprime_s(a_c,l_c,g_c,max_s(a_c,l_c,g_c));
            end
        end
    end

end %end age loop

% Remove the terminal value-function slice and map a' indices to grid values.
V(:,:,:,N_j+1) = [];
pol_aprime_val = a_grid(pol_aprime);

end %end function

function EV = f_expval(V_next,pi_l_vec,pi_g_vec)

n_a = size(V_next,1);
n_l = numel(pi_l_vec);
n_g = numel(pi_g_vec);

% Expected continuation value as a function of next-period assets a'.
EV = zeros(n_a,1);
for lprime=1:n_l
    for gprime=1:n_g
        EV = EV + V_next(:,lprime,gprime)*pi_l_vec(lprime)*pi_g_vec(gprime);
    end
end

end %end function

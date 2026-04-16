function [V,pol_s_idx,pol_s,pol_aprime_left_idx,pol_aprime_layer_idx,pol_aprime_val] = fun_vfi_given_ptype(n_a,n_s,n_l,n_g,N_j,...
    a_grid,s_grid,l_grid,g_grid,pi_l,pi_g,Params,opt)

V          = zeros(n_a,n_l,n_g,N_j+1);
pol_s_idx  = zeros(n_a,n_l,n_g,N_j);
pol_s      = zeros(n_a,n_l,n_g,N_j);
pol_aprime_left_idx = zeros(n_a,n_l,n_g,N_j);
pol_aprime_layer_idx = ones(n_a,n_l,n_g,N_j);

% Precompute asset grids in column/row form so that f_Return_cpu returns
% a matrix with rows indexed by a' and columns indexed by current assets a.
a_row      = a_grid';
if opt.gridinterplayer==1
    aprime_fine_grid = interp1((1:n_a)',a_grid,linspace(1,n_a,n_a+(n_a-1)*opt.ngridinterp)');
    aprime_col = aprime_fine_grid;
    %n_aprime = numel(aprime_fine_grid);
else
    aprime_fine_grid = a_grid;
    aprime_col = a_grid;
    %n_aprime = n_a;
end

% Backward induction over age.
for j_c = N_j:-1:1
    fprintf('Age = %d \n',j_c)
    V_next = V(:,:,:,j_c+1); % (a',l',g')
    V_s               = -Inf(n_a,n_l,n_g,n_s);
    pol_aprime_s_fine = zeros(n_a,n_l,n_g,n_s);

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
                EV = f_expval(V_next,pi_l_vec,pi_g_vec); 
                if opt.gridinterplayer==1
                    EVchoice = interp1(a_grid,EV,aprime_fine_grid); % (a',1), a' on finer grid
                else
                    EVchoice = EV; % (a',1), a' on coarse grid
                end
                % Current-period return matrix indexed by (a',a).
                RetMat_s_z = f_Return_cpu(s_val,aprime_col,a_row,l_val,g_val,...
                    Params.agej(j_c),Params.educ_i,Params.Jr,Params.r,Params.w,Params.ben,Params.pens,Params.gamma,Params.B_s);
                % Bellman RHS for the current (l,g,s) block.
                entireRHS = RetMat_s_z + Params.beta*EVchoice; % (a',a)
                [Vtemp,maxindex]   = max(entireRHS,[],1);  % (1,a)
                % Value and a' policy conditional on the current search choice.
                V_s(:,l_c,g_c,s_c) = Vtemp;
                pol_aprime_s_fine(:,l_c,g_c,s_c) = maxindex;
            end %end g
        end %end l
    end %end s

    % For the current age j, maximize over search effort.
    [V_jj,max_s] = max(V_s,[],4); % (a,l,g)
    V(:,:,:,j_c) = V_jj;
    pol_s_idx(:,:,:,j_c) = max_s; % index
    pol_s(:,:,:,j_c) = s_grid(max_s); % value

    % Recover the optimal a' index using the optimal search choice max_s.
    for a_c=1:n_a 
        for l_c=1:n_l
            for g_c=1:n_g
                fine_idx = pol_aprime_s_fine(a_c,l_c,g_c,max_s(a_c,l_c,g_c));
                if opt.gridinterplayer==1
                    % If interpolation is on, 
                    pol_aprime_left_idx(a_c,l_c,g_c,j_c) = floor((fine_idx-1)/(opt.ngridinterp+1))+1;
                    pol_aprime_layer_idx(a_c,l_c,g_c,j_c) = fine_idx-(opt.ngridinterp+1)*(pol_aprime_left_idx(a_c,l_c,g_c,j_c)-1);
                else
                    pol_aprime_left_idx(a_c,l_c,g_c,j_c) = fine_idx;
                end
            end
        end
    end

end %end age loop

% Remove the terminal value-function slice and map a' indices to grid values.
V(:,:,:,N_j+1) = [];
if opt.gridinterplayer==1
    pol_aprime_val = aprime_fine_grid((opt.ngridinterp+1)*(pol_aprime_left_idx-1)+pol_aprime_layer_idx);
else
    pol_aprime_val = a_grid(pol_aprime_left_idx);
end

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

function [V,pol_s,pol_aprime,StatDist,ValuesOnGrid,AllStats,AgeStats] = fun_solve2(Params,a_grid,s_grid,l_grid,g_grid,pi_l,pi_g,N_j,N_i)
% Solve the finite-horizon life-cycle problem without using the VFI Toolkit.
% ================================================
% State variables:
% a        Current endogenous state       a
% l        Employment state               semiz1
% g        Skill level                    semiz2
% educ_i   Permanent type                 i
% Choice variables:
% s        Search effort                  d
% aprime   Next-period endogenous state   aprime
% ================================================
% Outputs:
% V(a,l,g,j,i) stores the value function on the state grid.
% pol_s(a,l,g,j,i) stores the optimal search effort.
% pol_aprime(a,l,g,j,i) stores the optimal next-period asset choice in levels.
% StatDist(a,l,g,j,i) stores the agents distribution.
% ValuesOnGrid stores the evaluated objects on the state grid by ptype.
% AllStats stores unconditional means and means by permanent type.
% AgeStats stores conditional means by age and by age-permanent-type.

n_a = numel(a_grid);
n_s = numel(s_grid);
n_l = numel(l_grid);
n_g = numel(g_grid);

%% Value function iteration
disp('Start VFI in fun_solve2...')
V          = zeros(n_a,n_l,n_g,N_j,N_i);
pol_s_idx  = zeros(n_a,n_l,n_g,N_j,N_i);
pol_s      = zeros(n_a,n_l,n_g,N_j,N_i);
pol_aprime_idx = zeros(n_a,n_l,n_g,N_j,N_i);
pol_aprime = zeros(n_a,n_l,n_g,N_j,N_i);

Params_ii = Params;

% Loop over fixed types
for ii=1:N_i
    % pi_l(:,:,s_c,j_c,i_c)
    pi_l_ii = pi_l(:,:,:,:,ii);
    Params_ii.educ_i = Params.educ_i(ii);
    [V_ii,pol_s_idx_ii,pol_s_ii,pol_aprime_idx_ii,pol_aprime_ii] = fun_vfi_given_ptype(n_a,n_s,n_l,n_g,N_j,...
        a_grid,s_grid,l_grid,g_grid,pi_l_ii,pi_g,Params_ii);
    V(:,:,:,:,ii)          = V_ii;
    pol_s_idx(:,:,:,:,ii)  = pol_s_idx_ii;
    pol_s(:,:,:,:,ii)      = pol_s_ii;
    pol_aprime_idx(:,:,:,:,ii) = pol_aprime_idx_ii;
    pol_aprime(:,:,:,:,ii) = pol_aprime_ii;
end

%% Stationary distribution
disp('Start distribution in fun_solve2...')

StatDist = zeros(n_a,n_l,n_g,N_j,N_i);

% Initial condition
for ii=1:N_i
    % All start with zero assets, l_idx=1 (unemployed), mid skill
    StatDist(1,1,floor((n_g+1)/2),1,ii) = Params.omega_i(ii);
end

% Forward iterations
for ii=1:N_i
    for jj=1:1:N_j-1
        for g_c=1:n_g
            for l_c=1:n_l
                for a_c=1:n_a
                    s_c_opt = pol_s_idx(a_c,l_c,g_c,jj,ii);
                    aprime_c_opt = pol_aprime_idx(a_c,l_c,g_c,jj,ii);
                    prob_l = pi_l(l_c,:,s_c_opt,jj,ii);
                    prob_g = pi_g(g_c,:,l_c);
                    for lp_c=1:n_l
                        for gp_c=1:n_g
                            StatDist(aprime_c_opt,lp_c,gp_c,jj+1,ii)=...
                                StatDist(aprime_c_opt,lp_c,gp_c,jj+1,ii)+...
                                prob_l(lp_c)*prob_g(gp_c)*StatDist(a_c,l_c,g_c,jj,ii);
                        end
                    end
                end %end a
            end
        end
    end %end age
end %end ptype

% Convert cohort distributions into the cross-sectional distribution using
% the age-mass profile supplied in Params.mewj.
if ~isfield(Params,'mewj')
    error('fun_solve2 requires Params.mewj to be provided.')
end
if numel(Params.mewj)~=N_j
    error('Params.mewj must have length N_j.')
end
StatDist = StatDist.*reshape(Params.mewj,1,1,1,N_j,1);

%% Functions to evaluate
FnsToEvaluate.assets      = @(s,aprime,a,l,g) a;
FnsToEvaluate.assets_next = @(s,aprime,a,l,g) aprime;
FnsToEvaluate.search      = @(s,aprime,a,l,g) s;
FnsToEvaluate.empl        = @(s,aprime,a,l,g) l;
FnsToEvaluate.skill       = @(s,aprime,a,l,g) g;

%% Calculate values on the state grid
ValuesOnGrid = struct();
asset_grid_5d = repmat(reshape(a_grid,[n_a,1,1,1]),[1,n_l,n_g,N_j]);
empl_grid_5d = repmat(reshape(l_grid,[1,n_l,1,1]),[n_a,1,n_g,N_j]);
for ii=1:N_i
    name_i = sprintf('ptype%03d',ii);
    ValuesOnGrid.assets.(name_i) = asset_grid_5d;
    ValuesOnGrid.assets_next.(name_i) = pol_aprime(:,:,:,:,ii);
    ValuesOnGrid.search.(name_i) = pol_s(:,:,:,:,ii);
    ValuesOnGrid.empl.(name_i) = empl_grid_5d;
    ValuesOnGrid.skill.(name_i) = repmat(reshape(g_grid,[1,1,n_g,1]),[n_a,n_l,1,N_j]);
end

%% Calculate statistics and Life-cycle profiles

[AllStats,AgeStats] = fun_moments(StatDist,pol_s,pol_aprime,a_grid,l_grid,g_grid,FnsToEvaluate,N_j,N_i);

end %end function

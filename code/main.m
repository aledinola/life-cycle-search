clearvars,clc,close all

% Add the VFI Toolkit to the MATLAB path.
addpath(genpath('C:\Users\aledi\Documents\GitHub\VFIToolkit-matlab'))

%% Order of variables in the toolkit
% (d,a',a,semiz,z,j,i) = (s,a',a,l,g,age,educ)
% There is no exogenous state z in this model.

%% Demographics and state-space dimensions
% Model agents from age 20 to age 100, so there are 81 periods.

Params.agejshifter=19; % Age 20 minus one, useful for mapping model age to calendar age.
Params.J=100-Params.agejshifter; % =81, number of periods in the life cycle.
% Demographics
Params.agej=1:1:Params.J; % Model age index j = 1,2,...,J.
Params.Jr=46;             % Retirement starts at model age j = 46 (calendar age 65).

% Grid sizes to use
n_s     = 5;        % Search effort, decision d
n_a     = 501;      % Endogenous asset holdings
n_l     = 2;        % Employment states (l=0,1), semi-exogenous
n_g     = 3;       % Number of skills, semi-exogenous     
N_j     = Params.J; % Number of periods in finite horizon
N_i     = 2;        % Number of ptypes (education)

% semiz = (l,g) where l=employment, g=skill

%% Economic parameters

Params.w     = 1.0;   % Wage
Params.r     = 0.04;  % Real interest rate
Params.beta  = 0.988; % Discount factor
Params.csi   = 0.5;   % Search technology parameter
Params.gamma = 2.0;   % CRRA in utility consumption
Params.B_s   = 0.85;  % Search disutility parameter
Params.ben   = 0.1;   % Unemployment benefits
Params.pens  = 0.5;   % Pension
Params.educ_i  = [0.8,1.2]; % Fixed effect: low-high education
Params.sigma_i = [0.05,0.027]; % Prob of losing job: higher if low educ
Params.omega_i = [0.55,0.45]; % Fixed distribution of ptype
Params.mewj    = ones(1,Params.J)/N_j; % Cross-sectional mass by age

%% Search-effort choice

% Job-finding probability pi(s).
% Search-effort grid.
s_grid = linspace(0,1,n_s)';

%% Asset grid
a_min  = 0;
a_max  = 50;
a_grid = a_min+(a_max-a_min)*(linspace(0,1,n_a).^3)';

%% Skill grid and transition matrix
g_grid = [0.100, 0.500, 1.000]';

if numel(g_grid)~=n_g
    error('Length of g_grid is not correct.')
end

% When unemployed, skill can weakly depreciate by one rung.
pi_g_unemp = [1.000 0.0 0.0; 
              0.042 0.958 0.0;
              0.0   0.018 0.982];

% When employed, skill can weakly improve by one rung.
pi_g_emp   = [0.941 0.059 0.0;
              0.0   0.935 0.065;
              0.0   0.0 1.000];

% Transition matrix (g,g') depends on semi-exogenous state l=0,1
pi_g = zeros(n_g,n_g,n_l);
pi_g(:,:,1) = pi_g_unemp;
pi_g(:,:,2) = pi_g_emp;

% Basic probability checks for the skill transition matrices.
if any(abs(sum(pi_g_unemp,2)-1)>1e-12) || any(abs(sum(pi_g_emp,2)-1)>1e-12)
    error('Rows of pi_g must sum to one for every l')
end

%% Employment transitions

% Grid
l_grid = [0,1]'; % l=0 is unemployed, l=1 is employed
n_l    = numel(l_grid);

% Employment transitions depend on search effort through the job-finding
% probability, on permanent type, and on age. From retirement onward,
% employed agents separate with probability one, so employment is zero from
% age J_R+1 onward.
pi_l = zeros(n_l,n_l,n_s,N_j,N_i);
for i_c=1:N_i % ptype loop
    for j_c = 1:N_j % age loop
        if j_c >= Params.Jr
            sep_rate = 1.0;
        else
            sep_rate = Params.sigma_i(i_c);
        end
    for s_c = 1:n_s
        pi_l(:,:,s_c,j_c,i_c) = [1-f_pi_s(s_grid(s_c),Params.csi), f_pi_s(s_grid(s_c),Params.csi);
            sep_rate,                          1-sep_rate];
        if any(abs(sum(pi_l(:,:,s_c,j_c,i_c),2)-1)>1e-10)
            warning('Rows of pi_l do not sum to 1 for s=%d, age=%d and type=%d \n',s_c,j_c,i_c)
        end
    end
    end %end age loop
end %end ptype loop

% Basic probability checks for job-finding probabilities and transition rows.
if any(abs(sum(pi_l,2)-1)>1e-12,"all")
    error('Rows of pi_l must sum to one for every search choice.')
end

n_semiz = [n_l,n_g];
mid_g = floor((n_g+1)/2);
plot_j = min(5,N_j);
plot_i = 1;

%% Solve with toolkit-based implementation
start1 = tic;
[V, pol_s,pol_aprime,StatDist,ValuesOnGrid,AllStats,AgeStats] = fun_solve1(Params,a_grid,s_grid,l_grid,g_grid,pi_l,pi_g,N_j,N_i);
time1 = toc(start1);

% Check: pol_s(a,l,g,j,ptype) must be zero if l=1 (employed), for all
% a,g,j,ptype
if any(pol_s(:,2,:,:,:)>0,"all")
    warning('Optimal search effort >0 given employed, for some (a,g,j)')
end


%% Solve with direct MATLAB implementation
start2 = tic;
[V2,pol_s2,pol_aprime2,StatDist2,ValuesOnGrid2,AllStats2,AgeStats2] = fun_solve2(Params,a_grid,s_grid,l_grid,g_grid,pi_l,pi_g,N_j,N_i);
time2 = toc(start2);

fprintf('Time method 1 (GPU Toolkit): %f \n',time1)
fprintf('Time method 2 (GPU Toolkit): %f \n',time2)

err_V = max(abs(V(:)-V2(:)));
err_s = max(abs(pol_s(:)-pol_s2(:)));
err_aprime = max(abs(pol_aprime(:)-pol_aprime2(:)));
err_StatDist = max(abs(StatDist(:)-StatDist2(:)));
err_values_assets_ptype001 = max(abs(ValuesOnGrid.assets.ptype001(:)-ValuesOnGrid2.assets.ptype001(:)));
err_values_assets_ptype002 = max(abs(ValuesOnGrid.assets.ptype002(:)-ValuesOnGrid2.assets.ptype002(:)));
err_values_assets_next_ptype001 = max(abs(ValuesOnGrid.assets_next.ptype001(:)-ValuesOnGrid2.assets_next.ptype001(:)));
err_values_assets_next_ptype002 = max(abs(ValuesOnGrid.assets_next.ptype002(:)-ValuesOnGrid2.assets_next.ptype002(:)));
err_values_search_ptype001 = max(abs(ValuesOnGrid.search.ptype001(:)-ValuesOnGrid2.search.ptype001(:)));
err_values_search_ptype002 = max(abs(ValuesOnGrid.search.ptype002(:)-ValuesOnGrid2.search.ptype002(:)));
err_values_empl_ptype001 = max(abs(ValuesOnGrid.empl.ptype001(:)-ValuesOnGrid2.empl.ptype001(:)));
err_values_empl_ptype002 = max(abs(ValuesOnGrid.empl.ptype002(:)-ValuesOnGrid2.empl.ptype002(:)));
err_assets_mean = max(abs(AgeStats.assets.Mean(:)-AgeStats2.assets.Mean(:)));
err_assets_next_mean = max(abs(AgeStats.assets_next.Mean(:)-AgeStats2.assets_next.Mean(:)));
err_search_mean = max(abs(AgeStats.search.Mean(:)-AgeStats2.search.Mean(:)));
err_empl_mean = max(abs(AgeStats.empl.Mean(:)-AgeStats2.empl.Mean(:)));
err_assets_ptype001 = max(abs(AgeStats.assets.ptype001.Mean(:)-AgeStats2.assets.ptype001.Mean(:)));
err_assets_ptype002 = max(abs(AgeStats.assets.ptype002.Mean(:)-AgeStats2.assets.ptype002.Mean(:)));
err_search_ptype001 = max(abs(AgeStats.search.ptype001.Mean(:)-AgeStats2.search.ptype001.Mean(:)));
err_search_ptype002 = max(abs(AgeStats.search.ptype002.Mean(:)-AgeStats2.search.ptype002.Mean(:)));
err_empl_ptype001 = max(abs(AgeStats.empl.ptype001.Mean(:)-AgeStats2.empl.ptype001.Mean(:)));
err_empl_ptype002 = max(abs(AgeStats.empl.ptype002.Mean(:)-AgeStats2.empl.ptype002.Mean(:)));
err_allstats_assets = abs(AllStats.assets.Mean-AllStats2.assets.Mean);
err_allstats_assets_next = abs(AllStats.assets_next.Mean-AllStats2.assets_next.Mean);
err_allstats_search = abs(AllStats.search.Mean-AllStats2.search.Mean);
err_allstats_empl = abs(AllStats.empl.Mean-AllStats2.empl.Mean);
err_allstats_assets_ptype001 = abs(AllStats.assets.ptype001.Mean-AllStats2.assets.ptype001.Mean);
err_allstats_assets_ptype002 = abs(AllStats.assets.ptype002.Mean-AllStats2.assets.ptype002.Mean);
err_allstats_search_ptype001 = abs(AllStats.search.ptype001.Mean-AllStats2.search.ptype001.Mean);
err_allstats_search_ptype002 = abs(AllStats.search.ptype002.Mean-AllStats2.search.ptype002.Mean);
err_allstats_empl_ptype001 = abs(AllStats.empl.ptype001.Mean-AllStats2.empl.ptype001.Mean);
err_allstats_empl_ptype002 = abs(AllStats.empl.ptype002.Mean-AllStats2.empl.ptype002.Mean);
err_allstats_l0_assets = abs(AllStats.l0.assets.Mean-AllStats2.l0.assets.Mean);
err_allstats_l1_assets = abs(AllStats.l1.assets.Mean-AllStats2.l1.assets.Mean);
err_allstats_l0_skill = abs(AllStats.l0.skill.Mean-AllStats2.l0.skill.Mean);
err_allstats_l1_skill = abs(AllStats.l1.skill.Mean-AllStats2.l1.skill.Mean);
err_agestats_l0_assets = max(abs(AgeStats.l0.assets.Mean(:)-AgeStats2.l0.assets.Mean(:)));
err_agestats_l1_assets = max(abs(AgeStats.l1.assets.Mean(:)-AgeStats2.l1.assets.Mean(:)));
err_agestats_l0_skill = max(abs(AgeStats.l0.skill.Mean(:)-AgeStats2.l0.skill.Mean(:)));
err_agestats_l1_skill = max(abs(AgeStats.l1.skill.Mean(:)-AgeStats2.l1.skill.Mean(:)));

fprintf('err_V:        %f \n',err_V)
fprintf('err_s:        %f \n',err_s)
fprintf('err_aprime:   %f \n',err_aprime)
fprintf('err_StatDist: %f \n',err_StatDist)
fprintf('err_values_assets_ptype001:      %f \n',err_values_assets_ptype001)
fprintf('err_values_assets_ptype002:      %f \n',err_values_assets_ptype002)
fprintf('err_values_assets_next_ptype001: %f \n',err_values_assets_next_ptype001)
fprintf('err_values_assets_next_ptype002: %f \n',err_values_assets_next_ptype002)
fprintf('err_values_search_ptype001:      %f \n',err_values_search_ptype001)
fprintf('err_values_search_ptype002:      %f \n',err_values_search_ptype002)
fprintf('err_values_empl_ptype001:        %f \n',err_values_empl_ptype001)
fprintf('err_values_empl_ptype002:        %f \n',err_values_empl_ptype002)
fprintf('err_assets_mean:        %f \n',err_assets_mean)
fprintf('err_assets_next_mean:   %f \n',err_assets_next_mean)
fprintf('err_search_mean:        %f \n',err_search_mean)
fprintf('err_empl_mean:          %f \n',err_empl_mean)
fprintf('err_assets_ptype001:    %f \n',err_assets_ptype001)
fprintf('err_assets_ptype002:    %f \n',err_assets_ptype002)
fprintf('err_search_ptype001:    %f \n',err_search_ptype001)
fprintf('err_search_ptype002:    %f \n',err_search_ptype002)
fprintf('err_empl_ptype001:      %f \n',err_empl_ptype001)
fprintf('err_empl_ptype002:      %f \n',err_empl_ptype002)
fprintf('err_allstats_assets:          %f \n',err_allstats_assets)
fprintf('err_allstats_assets_next:     %f \n',err_allstats_assets_next)
fprintf('err_allstats_search:          %f \n',err_allstats_search)
fprintf('err_allstats_empl:            %f \n',err_allstats_empl)
fprintf('err_allstats_assets_ptype001: %f \n',err_allstats_assets_ptype001)
fprintf('err_allstats_assets_ptype002: %f \n',err_allstats_assets_ptype002)
fprintf('err_allstats_search_ptype001: %f \n',err_allstats_search_ptype001)
fprintf('err_allstats_search_ptype002: %f \n',err_allstats_search_ptype002)
fprintf('err_allstats_empl_ptype001:   %f \n',err_allstats_empl_ptype001)
fprintf('err_allstats_empl_ptype002:   %f \n',err_allstats_empl_ptype002)
fprintf('err_allstats_l0_assets:       %f \n',err_allstats_l0_assets)
fprintf('err_allstats_l1_assets:       %f \n',err_allstats_l1_assets)
fprintf('err_allstats_l0_skill:        %f \n',err_allstats_l0_skill)
fprintf('err_allstats_l1_skill:        %f \n',err_allstats_l1_skill)
fprintf('err_agestats_l0_assets:       %f \n',err_agestats_l0_assets)
fprintf('err_agestats_l1_assets:       %f \n',err_agestats_l1_assets)
fprintf('err_agestats_l0_skill:        %f \n',err_agestats_l0_skill)
fprintf('err_agestats_l1_skill:        %f \n',err_agestats_l1_skill)

if max(err_s,err_aprime)>1e-10    
    error('Mistake')
end

%% Figures
figure
plot(a_grid,a_grid,'--')
hold on
plot(a_grid,pol_aprime(:,1,mid_g,plot_j,plot_i))
hold on
plot(a_grid,pol_aprime(:,2,mid_g,plot_j,plot_i))
legend('45 line','Unemployed','Employed')
title('Policy aprime')

figure
plot(a_grid,a_grid,'--')
hold on
plot(a_grid,pol_aprime(:,2,1,plot_j,plot_i))
hold on
plot(a_grid,pol_aprime(:,2,n_semiz(2),plot_j,plot_i))
legend('45 line','Lowest skill','Highest skill')
title('Policy aprime')

figure
plot(a_grid,pol_s(:,1,mid_g,plot_j,plot_i))
hold on
plot(a_grid,pol_s(:,2,mid_g,plot_j,plot_i))
legend('Unemployed','Employed')
title('Policy search effort')


figure
plot(1:1:N_j,AgeStats.assets_next.ptype001.Mean)
hold on 
plot(1:1:N_j,AgeStats.assets_next.Mean)
hold on
plot(1:1:N_j,AgeStats.assets_next.ptype002.Mean)
title('Assets (a) by education type')
legend('Education: low type','Overall','Education: high type')

figure
plot(1:1:N_j,AgeStats.l0.assets_next.Mean)
hold on 
plot(1:1:N_j,AgeStats.assets_next.Mean)
hold on
plot(1:1:N_j,AgeStats.l1.assets_next.Mean)
title('Assets (a) by employment state')
legend('Unemployed','Overall','Employed')
 
figure
plot(1:1:N_j,AgeStats.search.ptype001.Mean)
hold on 
plot(1:1:N_j,AgeStats.search.Mean)
hold on
plot(1:1:N_j,AgeStats.search.ptype002.Mean)
title('Life cycle profile: Search effort (s)')
legend('Education: low type','Overall','Education: high type')

figure
plot(1:1:N_j,AgeStats.empl.ptype001.Mean)
hold on 
plot(1:1:N_j,AgeStats.empl.Mean)
hold on
plot(1:1:N_j,AgeStats.empl.ptype002.Mean)
title('Life cycle profile: Employment (l)')
legend('Education: low type','Overall','Education: high type')

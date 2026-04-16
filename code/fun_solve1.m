function [V,pol_s,pol_aprime,StatDist_out,ValuesOnGrid,AllStats,AgeStats,SimPanelValues,AgeStatsSim,ToolkitTimes] = fun_solve1(Params,a_grid,s_grid,l_grid,g_grid,pi_l,pi_g,N_j,N_i,opt)
% This function solves the lifecycle model using VFI toolkit
% VFI-Toolkit model: finite horizon semi-exo, no d1, no z
% The transition matrix of semiz is entered as an array in vfoptions.pi_semiz
% ================================================
% State variables:
% a        Current endogenous state       a
% l        Employment state               semiz1
% g        Skill level                    semiz2
% Choice variables:
% s        Search effort                  d
% aprime   Next-period endogenous state   aprime
% ================================================
% Outputs:
% V(a,l,g,j,i) stores the value function on the state grid.
% pol_s(a,l,g,j,i) stores the optimal search effort.
% pol_aprime(a,l,g,j,i) stores the optimal next-period asset choice in levels.
% StatDist_out(a,l,g,j,i) stores the agents distribution.
% Toolkit functions called by fun_solve1:
% - ValueFnIter_Case1_FHorz_PType
% - PolicyInd2Val_Case1_FHorz_PType
% - StationaryDist_Case1_FHorz_PType
% - ValuesOnGrid: EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1_PType
% - AllStats:     EvalFnOnAgentDist_AllStats_FHorz_Case1_PType
% - AgeStats:     LifeCycleProfiles_FHorz_Case1_PType

ToolkitTimes = struct();

n_a = length(a_grid);
n_s = length(s_grid);
n_l = length(l_grid);
n_g = length(g_grid);

%% Permanent types
PTypeDistName = {'omega_i'};
if ~isfield(Params,PTypeDistName{1})
    error('I did not find field %s in structure Params!',PTypeDistName{1})
end

if numel(Params.(PTypeDistName{1}))~=N_i
    error('No. of elements in PTypeDist is NOT equal to N_i')
end

%% Create big transition matrix pi_semiz(semiz,semiz'|d,i)
% where semiz = (l,g), d = search effort, i is permanent type

% pi_semiz is now (semiz,semiz'|d,j,i)
pi_semiz = create_pi_semiz(pi_l,pi_g,n_s,n_l,n_g,N_j,N_i);

% Switch to toolkit notation
n_d     = n_s;       % Num points for grid decision d
d_grid  = s_grid;    % Grid decision d --> semi-exo transitions
n_semiz = [n_l,n_g]; % semiz = (l,g)

% semiz_grid(semiz) OR semiz_grid(semiz,j) OR semiz_grid(semiz,j,ptype)
semiz_grid = [l_grid;g_grid]; % semiz = (l,g)
% I add age and type dimension to semiz (even if it is not needed)
semiz_grid_J = zeros(sum(n_semiz),N_j,N_i);
for ii=1:N_i
for jj=1:N_j
    semiz_grid_J(:,jj,ii) = semiz_grid;
end
end

vfoptions.n_semiz = n_semiz;
% Stacked grid for semiz:
vfoptions.semiz_grid      = semiz_grid_J; % semiz = (l,g)
vfoptions.pi_semiz        = pi_semiz; 
vfoptions.gridinterplayer = opt.gridinterplayer;
vfoptions.ngridinterp     = opt.ngridinterp;

simoptions.n_semiz    = vfoptions.n_semiz;
simoptions.semiz_grid = vfoptions.semiz_grid;
simoptions.pi_semiz   = vfoptions.pi_semiz;
simoptions.d_grid     = d_grid;
simoptions.gridinterplayer = vfoptions.gridinterplayer;
simoptions.ngridinterp     = vfoptions.ngridinterp;

n_z = 0;
z_grid = [];
pi_z   = [];

%% Now, create the return function 
Discs={'beta'};

% Use 'LifeCycleModel29_ReturnFn'
ReturnFn = @(s,aprime,a,l,g,agej,educ_i,Jr,r,w,ben,pens,gamma,B_s) ...
    f_Return(s,aprime,a,l,g,agej,educ_i,Jr,r,w,ben,pens,gamma,B_s);

%% Now solve the value function iteration problem, just to check that things are working before we go to General Equilbrium
 
disp('Start VFI in fun_solve1...')
vfoptions.verbose=1;
start_vfi = tic;
[V, Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,Discs,vfoptions);
ToolkitTimes.vfi = toc(start_vfi);

% V(a,semiz1,semiz2,age,ptype)
size(V.ptype001)
size(Policy.ptype001)

% Convert policy indexes to values
start_policy = tic;
PolicyVals=PolicyInd2Val_Case1_FHorz_PType(Policy,n_d,n_a,n_z,N_j,d_grid,a_grid,vfoptions);
ToolkitTimes.policy_values = toc(start_policy);

V_out      = zeros(n_a,n_semiz(1),n_semiz(2),N_j,N_i);
pol_s      = zeros(n_a,n_semiz(1),n_semiz(2),N_j,N_i);
pol_aprime = zeros(n_a,n_semiz(1),n_semiz(2),N_j,N_i);
PNames     = fieldnames(V);
for ii=1:N_i
    name_i = PNames{ii};
    V_out(:,:,:,:,ii)      = reshape(V.(name_i),[n_a,n_semiz(1),n_semiz(2),N_j]);
    pol_s(:,:,:,:,ii)      = reshape(PolicyVals.(name_i)(1,:,:,:),[n_a,n_semiz(1),n_semiz(2),N_j]);
    pol_aprime(:,:,:,:,ii) = reshape(PolicyVals.(name_i)(2,:,:,:),[n_a,n_semiz(1),n_semiz(2),N_j]);
end
V = V_out;

%% Distribution
disp('Start distribution in fun_solve1...')
% Initial distribution of agents at birth (j=1)
jequaloneDist = zeros(n_a,n_semiz(1),n_semiz(2)); 
jequaloneDist(1,1,floor((n_semiz(2)+1)/2)) = 1; 
AgeWeightsParamNames={'mewj'};

start_distribution = tic;
StatDist=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightsParamNames,PTypeDistName,Policy,n_d,n_a,n_z,N_j,N_i,pi_z,Params,simoptions);
ToolkitTimes.distribution = toc(start_distribution);

% Note: StatDist is the distribution conditional on ptype. To get the unconditional distribution, 
% we need to weight by the ptype distribution (Params.omega_i) 
PNames     = fieldnames(StatDist);
StatDist_out = zeros(n_a,n_semiz(1),n_semiz(2),N_j,N_i);
for ii=1:N_i
    name_i = PNames{ii};
    StatDist_out(:,:,:,:,ii) = Params.omega_i(ii)*reshape(StatDist.(name_i),[n_a,n_semiz(1),n_semiz(2),N_j]);
end

%% Functions to evaluate and conditional restrictions
FnsToEvaluate.assets=@(s,aprime,a,l,g) a; % a: current assets
FnsToEvaluate.assets_next=@(s,aprime,a,l,g) aprime; % a': next-period assets
FnsToEvaluate.search=@(s,aprime,a,l,g) s; % s: search effort
FnsToEvaluate.empl  =@(s,aprime,a,l,g) l; % l: employment (0-1)
FnsToEvaluate.skill =@(s,aprime,a,l,g) g; % g: skill

simoptions.conditionalrestrictions.l0 = @(s,aprime,a,l,g) l==0; % unemployed
simoptions.conditionalrestrictions.l1 = @(s,aprime,a,l,g) l==1; % employed

%% Calculate values on the state grid
start_values_on_grid = tic;
ValuesOnGrid = EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1_PType(Policy,FnsToEvaluate,...
    Params,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid,simoptions);
ToolkitTimes.values_on_grid = toc(start_values_on_grid);

%% Calculate statistics unconditional on age
% NOTE: We need only averages, so we add option to simoptions
simoptions.whichstats = zeros(7,1);
simoptions.whichstats(1) = 1; % only Mean

start_all_stats = tic;
AllStats = EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StatDist,Policy,FnsToEvaluate,...
    Params,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid,simoptions);
ToolkitTimes.all_stats = toc(start_all_stats);

%% Calculate the life-cycle profiles using the distribution StatDist
start_age_stats = tic;
AgeStats=LifeCycleProfiles_FHorz_Case1_PType(StatDist,Policy,FnsToEvaluate,...
    Params,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid,simoptions);
ToolkitTimes.age_stats = toc(start_age_stats);

%% Generate a simulated panel and from there life cycle profiles

simoptions_panel = simoptions;
simoptions_panel.semiz_grid = vfoptions.semiz_grid;
simoptions_panel.pi_semiz = vfoptions.pi_semiz;
simoptions_panel.alreadygridvals_semiexo = 0;
simoptions_panel.numbersims = 5000;
simoptions_panel.simperiods = N_j;
start_panel = tic;
SimPanelValues=SimPanelValues_FHorz_Case1_PType(jequaloneDist,PTypeDistName,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid,pi_z,simoptions_panel);
AgeStatsSim = panel_age_stats_ptype(SimPanelValues,Params.(PTypeDistName{1}),N_i);
ToolkitTimes.sim_panel_and_age_stats = toc(start_panel);
ToolkitTimes.total_subparts = ToolkitTimes.vfi + ToolkitTimes.policy_values + ToolkitTimes.distribution + ...
    ToolkitTimes.values_on_grid + ToolkitTimes.all_stats + ToolkitTimes.age_stats + ToolkitTimes.sim_panel_and_age_stats;

end %end function

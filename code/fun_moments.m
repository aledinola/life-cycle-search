function [AllStats,AgeStats] = fun_moments(StatDist,pol_s,pol_aprime,a_grid,l_grid,g_grid,FnsToEvaluate,N_j,N_i)
% Compute means of user-supplied functions from the cross-sectional distribution.
%
% Inputs
% StatDist(a,l,g,j,i) is the cross-sectional distribution over states and
% permanent types.
% pol_s(a,l,g,j,i) and pol_aprime(a,l,g,j,i) are the policy functions.
% FnsToEvaluate is a structure of function handles. Each field is a
% function of the form fn_handle(s,aprime,a,l,g), evaluated on the full
% state grid.
%
% Outputs
% AllStats is a structure with one field per function in FnsToEvaluate.
% For each function name fn_name:
% AllStats.(fn_name).Mean is the unconditional cross-sectional mean.
% AllStats.(fn_name).ptype001.Mean, AllStats.(fn_name).ptype002.Mean, ...
% are the means conditional on permanent type.
% AllStats.l0.(fn_name).Mean and AllStats.l1.(fn_name).Mean are the means
% conditional on employment state, with analogous by-ptype fields nested
% underneath each condition.
%
% AgeStats is a structure with one field per function in FnsToEvaluate.
% For each function name fn_name:
% AgeStats.(fn_name).Mean is a 1-by-N_j vector of means conditional on age.
% AgeStats.(fn_name).ptype001.Mean, AgeStats.(fn_name).ptype002.Mean, ...
% are 1-by-N_j vectors of means conditional on age and permanent type.
% AgeStats.l0.(fn_name).Mean and AgeStats.l1.(fn_name).Mean are 1-by-N_j
% vectors conditional on age and employment state, again with analogous
% by-ptype fields nested underneath each condition.

%% Calculate values on state grid for functions to evaluate
fn_names = fieldnames(FnsToEvaluate);
ptype_names = local_ptype_names(N_i);
cond_names = {'l0','l1'};

n_a = numel(a_grid);
n_l = numel(l_grid);
n_g = numel(g_grid);

[A,L,G] = ndgrid(a_grid,l_grid,g_grid);
ValuesOnGrid = struct();
for ff=1:numel(fn_names)
    fn_name = fn_names{ff};
    fn_handle = FnsToEvaluate.(fn_name);
    ValuesOnGrid.(fn_name) = zeros(n_a,n_l,n_g,N_j,N_i);
    for ii=1:N_i
        for jj=1:N_j
            ValuesOnGrid.(fn_name)(:,:,:,jj,ii) = fn_handle( ...
                pol_s(:,:,:,jj,ii),pol_aprime(:,:,:,jj,ii),A,L,G);
        end
    end
end

%% Conditional restriction masks on employment
L_full = reshape(L,[n_a,n_l,n_g,1,1]);
CondMasks = struct();
CondMasks.l0 = repmat(L_full==0,[1,1,1,N_j,N_i]);
CondMasks.l1 = repmat(L_full==1,[1,1,1,N_j,N_i]);

%% Calculate the mean unconditional and conditional on ptype, for all functions
AllStats = struct();
for ff=1:numel(fn_names)
    fn_name = fn_names{ff};
    vals = ValuesOnGrid.(fn_name);
    total_mass = sum(StatDist,'all');
    if total_mass>0
        AllStats.(fn_name).Mean = sum(vals.*StatDist,'all')/total_mass;
    else
        AllStats.(fn_name).Mean = NaN;
    end

    for ii=1:N_i
        ptype_mass = sum(StatDist(:,:,:,:,ii),'all');
        if ptype_mass>0
            AllStats.(fn_name).(ptype_names{ii}).Mean = ...
                sum(vals(:,:,:,:,ii).*StatDist(:,:,:,:,ii),'all')/ptype_mass;
        else
            AllStats.(fn_name).(ptype_names{ii}).Mean = NaN;
        end
    end

    for cc=1:numel(cond_names)
        cond_name = cond_names{cc};
        cond_dist = StatDist.*CondMasks.(cond_name);
        cond_weights_group = zeros(1,N_i);
        cond_means_group = NaN(1,N_i);

        for ii=1:N_i
            ptype_mass = sum(StatDist(:,:,:,:,ii),'all');
            cond_ptype_dist = cond_dist(:,:,:,:,ii);
            cond_ptype_mass = sum(cond_ptype_dist,'all');
            if cond_ptype_mass>0
                AllStats.(cond_name).(fn_name).(ptype_names{ii}).Mean = ...
                    sum(vals(:,:,:,:,ii).*cond_ptype_dist,'all')/cond_ptype_mass;
                cond_means_group(ii) = AllStats.(cond_name).(fn_name).(ptype_names{ii}).Mean;
            else
                AllStats.(cond_name).(fn_name).(ptype_names{ii}).Mean = NaN;
            end
            if ptype_mass>0
                cond_weights_group(ii) = cond_ptype_mass/ptype_mass;
            end
        end

        if sum(cond_weights_group)>0
            AllStats.(cond_name).(fn_name).Mean = ...
                sum(cond_weights_group.*cond_means_group,'omitnan')/sum(cond_weights_group);
        else
            AllStats.(cond_name).(fn_name).Mean = NaN;
        end
    end
end

%% Calculate the mean by age, for all functions
AgeStats = struct();
for ff=1:numel(fn_names)
    fn_name = fn_names{ff};
    AgeStats.(fn_name).Mean = zeros(1,N_j);
    for ii=1:N_i
        AgeStats.(fn_name).(ptype_names{ii}).Mean = zeros(1,N_j);
    end
    for cc=1:numel(cond_names)
        cond_name = cond_names{cc};
        AgeStats.(cond_name).(fn_name).Mean = zeros(1,N_j);
        for ii=1:N_i
            AgeStats.(cond_name).(fn_name).(ptype_names{ii}).Mean = zeros(1,N_j);
        end
    end
end

for ff=1:numel(fn_names)
    fn_name = fn_names{ff};
    vals = ValuesOnGrid.(fn_name);

    for jj=1:N_j
        age_mass = sum(StatDist(:,:,:,jj,:),'all');
        if age_mass>0
            AgeStats.(fn_name).Mean(jj) = sum(vals(:,:,:,jj,:).*StatDist(:,:,:,jj,:),'all')/age_mass;
        else
            AgeStats.(fn_name).Mean(jj) = NaN;
        end
    end

    %% Calculate the mean by (age-ptype), for all functions
    for ii=1:N_i
        for jj=1:N_j
            ptype_mass = sum(StatDist(:,:,:,jj,ii),'all');
            if ptype_mass>0
                AgeStats.(fn_name).(ptype_names{ii}).Mean(jj) = ...
                    sum(vals(:,:,:,jj,ii).*StatDist(:,:,:,jj,ii),'all')/ptype_mass;
            else
                AgeStats.(fn_name).(ptype_names{ii}).Mean(jj) = NaN;
            end
        end
    end

    %% Calculate the mean by employment condition, and by (age-condition-ptype)
    for cc=1:numel(cond_names)
        cond_name = cond_names{cc};
        cond_dist = StatDist.*CondMasks.(cond_name);
        for ii=1:N_i
            for jj=1:N_j
                cond_ptype_age_dist = cond_dist(:,:,:,jj,ii);
                cond_ptype_age_mass = sum(cond_ptype_age_dist,'all');
                if cond_ptype_age_mass>0
                    AgeStats.(cond_name).(fn_name).(ptype_names{ii}).Mean(jj) = ...
                        sum(vals(:,:,:,jj,ii).*cond_ptype_age_dist,'all')/cond_ptype_age_mass;
                else
                    AgeStats.(cond_name).(fn_name).(ptype_names{ii}).Mean(jj) = NaN;
                end
            end
        end

        for jj=1:N_j
            cond_weights_group = zeros(1,N_i);
            cond_means_group = NaN(1,N_i);
            for ii=1:N_i
                ptype_mass = sum(StatDist(:,:,:,:,ii),'all');
                cond_ptype_age_dist = cond_dist(:,:,:,jj,ii);
                cond_ptype_age_mass = sum(cond_ptype_age_dist,'all');
                cond_means_group(ii) = AgeStats.(cond_name).(fn_name).(ptype_names{ii}).Mean(jj);
                if ptype_mass>0
                    cond_weights_group(ii) = cond_ptype_age_mass/ptype_mass;
                end
            end
            if sum(cond_weights_group)>0
                AgeStats.(cond_name).(fn_name).Mean(jj) = ...
                    sum(cond_weights_group.*cond_means_group,'omitnan')/sum(cond_weights_group);
            else
                AgeStats.(cond_name).(fn_name).Mean(jj) = NaN;
            end
        end
    end
end

end %end function

function ptype_names = local_ptype_names(N_i)

ptype_names = cell(1,N_i);
for ii=1:N_i
    if ii<10
        ptype_names{ii} = ['ptype00',num2str(ii)];
    elseif ii<100
        ptype_names{ii} = ['ptype0',num2str(ii)];
    else
        ptype_names{ii} = ['ptype',num2str(ii)];
    end
end

end %end function

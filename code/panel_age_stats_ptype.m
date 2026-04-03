function AgeStatsSim = panel_age_stats_ptype(SimPanelValues,PTypeDist,N_i)
% Compute age profiles from a simulated panel with fixed-type blocks.
% Each field of SimPanelValues is [N_j,numbersims], with columns ordered by
% permanent type in the same block order used by SimPanelValues_FHorz_Case1_PType.

fn_names = fieldnames(SimPanelValues);
ptype_names = local_ptype_names(N_i);
ptype_counts = local_ptype_counts(PTypeDist,size(SimPanelValues.(fn_names{1}),2));

AgeStatsSim = struct();
for ff=1:numel(fn_names)
    fn_name = fn_names{ff};
    values = SimPanelValues.(fn_name);
    AgeStatsSim.(fn_name).Mean = mean(values,2,'omitnan')';

    start_idx = 1;
    for ii=1:N_i
        stop_idx = start_idx + ptype_counts(ii) - 1;
        if ptype_counts(ii)>0
            AgeStatsSim.(fn_name).(ptype_names{ii}).Mean = mean(values(:,start_idx:stop_idx),2,'omitnan')';
        else
            AgeStatsSim.(fn_name).(ptype_names{ii}).Mean = NaN(1,size(values,1));
        end
        start_idx = stop_idx + 1;
    end
end

end %end function

function ptype_counts = local_ptype_counts(PTypeDist,numbersims)

ptype_counts = floor(PTypeDist(:)'*numbersims);
extra_sims = numbersims-sum(ptype_counts);
ptype_counts(1:extra_sims) = ptype_counts(1:extra_sims) + 1;

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

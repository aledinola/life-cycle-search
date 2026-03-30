function prob=f_SemiExoStateFn(l,g,lprime,gprime,s, csi,sigma,p0,p1)
% Variable Description                    Toolkit
% ================================================
% s        Search effort                  d
% l        Employment state 0,1           semiz1
% g        Skill level                    semiz2
% ================================================
% Transitions
% Prob(l,l') depends on search effort
% Prob(g,g') depends on l and on age

prob=-1; % Just a placeholder (one that will cause errors if not overwritten)
prob_l=-1;
prob_g=Inf;

%% Transition from l to l'
if l==0
    if lprime==0
        % l=0 and l'=0
        prob_l = 1-pi_s(s,csi);
    elseif lprime==1
        % l=0 and l'=1: job finding probability
        prob_l = pi_s(s,csi);
    end
elseif l==1
    if lprime==0
        % l=1 and l'=0: job loss probability
        prob_l = sigma;
    elseif lprime==1
        % l=1 and l'=1
        prob_l = 1-sigma;
    end
end %end l

%% Transition from g to g'


%% Combine the l-probability with the g-probability
prob=prob_l*prob_g;



end  %end function

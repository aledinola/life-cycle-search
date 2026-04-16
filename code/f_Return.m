function F = f_Return(s,aprime,a,l,g,z3, agej,educ_i,Jr,r,w,ben,pens,gamma,B_s)
% Period return for the toolkit-based solver.
% Inputs follow the toolkit convention: choices (s,aprime) and current 
% states (a,l,g) are compulsory. 
% ================================================
% Variable Description                    Toolkit name
% ================================================
% s        Search effort                  d
% aprime   Next-period endogenous state   aprime
% a        Current endogenous state       a
% l        Employment state               semiz1
% g        Skill level                    semiz2
% z3       Check                          semiz3
% ================================================

if agej<Jr % Working age
    c = (1+r)*a + educ_i*(w*g*(l==1) + ben*(l==0)) - aprime + 0*z3;
else % Retirement
    c = (1+r)*a + educ_i*pens - aprime;
end

F=-Inf;
if c>0
    if gamma==1
        F = log(c) - B_s*s;
    else
        F = c^(1-gamma)/(1-gamma) - B_s*s;
    end
end

end %end function

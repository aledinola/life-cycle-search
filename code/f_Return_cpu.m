function F = f_Return_cpu(s,aprime,a,l,g, agej,educ_i,Jr,r,w,ben,pens,gamma,B_s)
% Period return for the direct MATLAB solver.
% This function mirrors f_Return.m, but it accepts array inputs so it can
% evaluate the full (a',a) return matrix in one call.
% ================================================
% s        Search effort                  d
% aprime   Next-period endogenous state   aprime
% a        Current endogenous state       a
% l        Employment state               semiz1
% g        Skill level                    semiz2
% ================================================

if agej<Jr % Working age
    c = (1+r)*a + educ_i*(w*g*(l==1) + ben*(l==0)) - aprime;
else % Retirement
    c = (1+r)*a + educ_i*pens - aprime;
end

F=-Inf(size(c));
c_pos = c>0;

if gamma==1
    F(c_pos) = log(c(c_pos)) - B_s*s;
else
    F(c_pos) = c(c_pos).^(1-gamma)/(1-gamma) - B_s*s;
end

end %end function

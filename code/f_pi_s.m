function pi_s = f_pi_s(s,csi)
% f_pi_s gives job finding probability as function of search effort s.

pi_s = min(csi.*s,1);

end
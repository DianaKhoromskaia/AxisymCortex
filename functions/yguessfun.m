function y = yguessfun(s,sol,sfun)
% This function generates an initial guess for the bvp4c-solver from the
% old solution sol mapped onto the updated shape via sfun(s) (old arc length as a function
% of new arc length). 

y = deval(sol, sfun(s));
end


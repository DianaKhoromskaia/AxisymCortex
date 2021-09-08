function y = yguessfun_initial(s, sol, factor)
% This function generates an initial guess for the bpv4c-solver from sol, scaled by factor.
 
y = factor*deval(sol, s);
end


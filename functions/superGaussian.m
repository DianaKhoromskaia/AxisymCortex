function f = superGaussian(x, x0, A, Aconst, s, P)
% super-Gaussian function with power P

f = A*exp(-(((x-x0*ones(1,length(x))).*(x-x0*ones(1,length(x))))/(2*s^2)).^P)+Aconst;

end


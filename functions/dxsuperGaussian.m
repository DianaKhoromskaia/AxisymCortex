function f = dxsuperGaussian(x, x0, A, s, P)
% First derivative of super-Gaussian with power P

f =-(A*(2^(1 - P))*P/(s^(2*P)))*exp(-2^(-P)*((x-x0).*(x-x0)/(s^2)).^P).*(x-x0).^(2*P-1);

end


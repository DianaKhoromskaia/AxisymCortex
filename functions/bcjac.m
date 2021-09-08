function [dbcdya, dbcdyb, dbcdp] = bcjac(ya, yb, varpar,  C1, C2, C, dsC,Psi, X, Z, X0, L, zeta, dszeta, eta, etab, etac, etacb, eps1, xintegral, fext, kappa, C0, delt, xi)
% Jacobian of the residual function for bvp4c
% size of matrix = (nBC, nreg x nequations)

dbcdya = zeros(12,11); %Jacobian wrt boundary values at s=0
dbcdyb = zeros(12,11); %Jacobian wrt boundary values at s=L
dbcdp = zeros(12,1); %Jacobian wrt free parameters

%partial derivatives wrt ya
dbcdya(1,1) = 1;
dbcdya(2,4) = 1;
dbcdya(3,6) = 1;
dbcdya(4,7) = 1;
dbcdya(5,8) = 1;
dbcdya(6,9) = 1;
dbcdya(7,10) = 1;

dbcdya(12,11) = 1;

%partial derivatives wrt yb
dbcdyb(8,1) = 1;
dbcdyb(9,4) = 1;
dbcdyb(10,7) = 1;
dbcdyb(11,6) = 1;

end


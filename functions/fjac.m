function [dfdv,dfdvar] = fjac(s, v, varpar, C1, C2, C, dsC, Psi, X, Z, X0, L, zeta, dszeta, eta, etab, etac, etacb, eps1, xintegral, fext, kappa, C0, delt, xi)
%Jacobian of RHS for the bvp4c solver

dfdv = zeros(11,11); %Jacobian wrt vector of unknowns
dfdvar = zeros(11,1); %Jacobian wrt free parameters

switch s
    case 0. % for s=0 (south pole)
        %partial derivatives wrt v-components
        dfdv(1,2) = 1;
        
        dfdv(4,3) = -0.5*(C1(0.).^2 + C2(0.).^2);
        dfdv(4,5) = -0.5/etacb;
        
        dfdv(6,2) = etab*C(0.);
        dfdv(6,3) = 0.5*etab*(C(0.).^2);
        dfdv(6,5) = -0.5*(C1(0.).^2 + C2(0.).^2);
        
        dfdv(9,2) = 1;
        dfdv(9,3) = C2(0);
        
        dfdv(11,3) = C2(0);
        
        %partial derivatives wrt to parameters
        dfdvar(6,1) = -0.5;
    
    otherwise % for s \in (0,L]
        
        %partial derivatives wrt v-components
        dfdv(1,2) = 1;
        
        dfdv(2,1) = (cos(Psi(s))./X(s)).^2 - ((eta-etab)/(eta+etab))*C1(s).*C2(s);
        dfdv(2,2) = -cos(Psi(s))./X(s);
        dfdv(2,3) = -dsC(s);
        dfdv(2,4) = -(etab*C(s) + eta*(C2(s)-C1(s)))/(eta+etab);
        dfdv(2,5) = dsC(s)/(eta+etab);
        
        dfdv(3,4) = 1;
        
        dfdv(4,1) = dsC(s);
        dfdv(4,3) = -(C1(s).^2 + C2(s).^2);
        dfdv(4,4) = -cos(Psi(s))./X(s);
        dfdv(4,5) = -1/etacb;
        
        dfdv(5,6) = 1;
        
        dfdv(6,1) = (eta*(C1(s)-C2(s)) + etab*C(s)).*cos(Psi(s))./X(s) - C1(s).*(etab-eta).*(cos(Psi(s)))./X(s);
        dfdv(6,2) = eta*(C2(s)-C1(s)) + etab*C(s) -C1(s).*(eta+etab);
        dfdv(6,3) = eta*(C1(s)-C2(s)).*(C1(s)-C2(s)) + etab*(C(s).*C(s)) -C1(s).*(2*eta*C2(s) +(etab-eta).*C(s));
        dfdv(6,5) = - (C1(s).^2 + C2(s).^2) + C1(s).*C2(s);
        dfdv(6,10) = -1/(X(s).*X(s));
        
        dfdv(7,3) = 2*pi*X(s);
        
        dfdv(8,3) = X(s).*(C(s).*(Z(s)-X0) - cos(Psi(s)))/xintegral;
        
        dfdv(9,2) = 1;
        dfdv(9,3) = C2(s);
        
        dfdv(11,3) = C2(s);
        
        %partial derivatives wrt to parameters
        dfdvar(6,1) = -0.5;

end
end


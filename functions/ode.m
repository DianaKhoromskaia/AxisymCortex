function dvds = ode(s, v, varpar, C1, C2, C, dsC, Psi, X, Z, X0, L, zeta, dszeta, eta, etab, etac, etacb, eps1, xintegral, fext, kappa, C0, delt, xi)
% dsvs is the RHS of the ode describing the force balance of an active viscous cortex
% for a given shape at time t, and auxiliary equations required for the
% dynamics.
% The vector of unknowns is
% v = (v_s, \partial_s v_s, v_n, \partial_s v_n, \bar{m}_s^s, t_n^s, dV(s), dX(s), d((snew-sold)_{Lagrange})/dt, I1(s), d((snew-sold)_{Euler})/dt)
% The components of dvds are:
%      - dvds(1:2): tangential force balance (D.3), written as a second order system for v_s
%      - dvds(3:4): constitutive equation (D.1), written as a second order system for v_n
%      - dvds(5:6): normal force balance (D.4), written as a second order
%      system for \bar{m}_s^s. Here, we have rewritten the NFB using the
%      expression (see Ref. [72])
%      	t_n ^s = t^s_s \tan\psi - \frac{1}{2}\frac{x}{\cos\psi}P +
%      	\frac{1}{x \cos\psi} I1,
%      which contains the partially integrated external force I1.
%      - dvds(7): equation for partial rate of volume change dV; with the
%      two corresponding b.c.'s this ensures incompressibility via
%      adjustment of the pressure P
%      - dvds(8): equation for partial centroid speed dX; in general, the free parameter fc (which enters the equations as a constant external force f^{ext}=fc \mathbf{e}_z) 
%      allows to fix the centroid position in time; with up-down symmetry this is not required, hence fc is not a parameter here and set to fc=0. 
%      - dvds(9): equation for the rate of change of the arc length parameter
%      on the surface at time t and the surface displaced by (v_s,v_n) (Lagrangian approach)
%      - dvds(10): equation for the partial integrated external force
%      I1(s); here I1(s)=0 for all s. 
%      - dvds(11): equation for the rate of change of the arc length parameter
%      on the surface at time t and the surface displaced by (0,v_n) (Eulerian approach)
% The solution is obtained on the interval s \in [0,L], where s=L is the equator of the shape, 
% with RHS specified separately for the south pole s=0 as the analytical limit of the expressions. 

P = varpar(1);
fc = 0;

ds2vs = -cos(Psi(s)).*(v(2,:)-cos(Psi(s)).*v(1,:)./X(s))./X(s) - (eta-etab)*C1(s).*C2(s).*v(1,:)/(eta+etab) - dsC(s).*v(3,:) - (etab*C(s)+eta*(C2(s)-C1(s))).*v(4,:)/(eta+etab) - dszeta(s)/(eta+etab) - sin(Psi(s)).*(fext(s)+fc)./(eta+etab) + dsC(s).*v(5,:)./(eta+etab);
ds2vn = -cos(Psi(s)).*v(4,:)./X(s) - v(3,:).*(C1(s).^2 + C2(s).^2) + v(1,:).*dsC(s) - v(5,:)/etacb;
Cijtij = v(1,:).*(cos(Psi(s))./X(s)).*(eta*(C1(s)-C2(s)) + etab*C(s)) + v(2,:).*(eta*(C2(s)-C1(s)) + etab*C(s)) + v(3,:).*(eta*(C1(s)-C2(s)).*(C1(s)-C2(s)) + etab*C(s).*C(s)) + C(s).*zeta(s);
vkk = v(2,:) + v(1,:).*cos(Psi(s))./X(s) + C(s).*v(3,:);
tss = 2*eta*(v(2,:)+C2(s).*v(3,:)) + (etab-eta)*vkk + zeta(s) - v(5,:).*C2(s); 
ds2mss = - P/2 + Cijtij - C1(s).*tss - v(5,:).*(C1(s).^2 + C2(s).^2);

% RHS for s \in (0,L]:
dvds = [v(2,:);
    ds2vs;
    v(4,:);
    ds2vn;
    v(6,:);
    ds2mss - v(10,:)./(X(s).*X(s)) + cos(Psi(s)).*(fext(s)+fc);
    2*pi*X(s).*v(3,:);
    X(s).*v(3,:).*(C(s).*(Z(s)-X0) - cos(Psi(s)))/xintegral;
    (v(2,:) + C2(s).*v(3,:));
    X(s).*(fext(s)+fc);
    (C2(s).*v(3,:))];

% RHS for s=0 (south pole):
indices = find(~s);
zervec = zeros(size(s(indices)));
onevec = ones(size(s(indices)));

dvds(:,indices) =   [v(2,indices);
    zervec;
    zervec;
    0.5*(-v(3,indices).*(C1(0).^2 + C2(0).^2) - v(5,indices)/etacb);
    zervec;
    0.5*(- P*onevec + (fext(0)+fc)*onevec + C(0).*(etab*(2*v(2,indices) + C(0).*v(3,indices)) + zeta(0.)*onevec) - v(5,indices).*(C1(0).^2 + C2(0).^2));
    zervec;
    zervec;
    (v(2,indices) + C2(0)*v(3,indices));
    zervec;
    C2(0)*v(3,indices)];
    
end


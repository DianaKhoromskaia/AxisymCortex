function [C1new, C2new, Cnew, Xnew, Psinew, Znew, dsC1new, dsC2new, dsCnew, zetanew, dszetanew, xintegral, x0, s0new] = evolvefunctions(svec, vs, dsvs, vn, dsvn, mss, tns, ds2vn, C1, C2, C, dsC1, dsC2, dsC, X, Psi, Z, snewfun, snewvec, snewfunEu, snewvecEu, t, dt, L, Lnew, eps1, eta, etac, etacb, gamma, x0fac, A, sigma0, pow, s0, zeta, npoints)
% This function evolves the shape in the Lagrangian approach, where material points move according to the full flow field
% \mathbf{v}=v^s \mathbf{e}_s + v^n \mathbf{n} (see Eqn. (D.9)); 
% The updated shape descriptors are saved as spline interpolants over the arc length parameter snew. 
% The tension profile moves with the normal component v_n \mathbf{n} only, in the Eulerian approach to match the 3D simulations.  
    
    % Define grids on the displaced surface:
    svec3 = [0. svec((svec > eps1))]; % adjust solver grid: remove grid points which are too close to the pole
    svec2 = svec3(2:end-1); % vector of svec3 without endpoints
    svec21 = svec3(2:end); % vector of svec3 without SP
    snew = [0. snewfun(svec2) Lnew]; % adjusted solver grid mapped onto new arc length via Lagrange update
    snewEu = [0. snewfunEu(svec2) Lnew]; % adjusted solver grid mapped onto new arc length via Euler update

    % Update the shape:
    % The expression below follow from the time derivative of (D.5) and of the normal vector \mathbf{n}  
    Xnew = griddedInterpolant( snew, [0. (X(svec21) + dt*(sin(Psi(svec21)).*vn(svec21) + cos(Psi(svec21)).*vs(svec21)))] , 'spline');
    Psinew = griddedInterpolant( snew, [0. (Psi(svec2) + dt*( -dsvn(svec2) + C2(svec2).*vs(svec2))) pi/2], 'spline');
    Znew = griddedInterpolant(snew, Z(svec3)+dt*(-cos(Psi(svec3)).*vn(svec3) + sin(Psi(svec3)).*vs(svec3)), 'spline');
    
    % Update the curvatures and curvature derivatives:
    
    % C=C_k^k and its derivative are updated using the constitutive equation (D.1);
    % The Jacobian factor is jac=\partial_s snew and is used for the
    % curvature derivative, when expressed in the new arc length.
    Cnew0 = C(0) + dt*mss(0.)/etacb;
    CnewL = C(L) + dt*mss(L)/etacb;
    jac = abs(ones(size(svec2)) + dt*(dsvs(svec2) + vn(svec2).*C2(svec2))); 
    Cnew = griddedInterpolant( snew, [Cnew0 (C(svec2) + dt*( mss(svec2)/etacb) ) CnewL], 'spline');
    dsCnew = griddedInterpolant( snew, [0. (dsC(svec2) + dt*( tns(svec2)/etacb))./jac 0.], 'spline');
    
    % C1=C^{\phi}^{\phi} is updated using the corotational time derivative
    % of the curvature tensor (see Ref. [40]); its spatial derivative is
    % obtained from the relation \partial_s C^{\phi}^{\phi} =
    % (C_s^s-C^{\phi}^{\phi})*cos(psi)/x, which holds for axisymmetric
    % surfaces. C2=C^s_s is obtained from C2=C-C1.
    C1new = griddedInterpolant(snew, [Cnew0/2 C1(svec21)+dt*(-dsvn(svec21).*cos(Psi(svec21))./X(svec21) -vn(svec21).*C1(svec21).*C1(svec21) + vs(svec21).*dsC1(svec21) )], 'spline');  
    C2new = griddedInterpolant(snew, Cnew(snewfun(svec3))-C1new(snewfun(svec3)), 'spline');
    
    dsC1new = griddedInterpolant( snew, [0. cos(Psinew(snewfun(svec2))).*(C2new(snewfun(svec2))-C1new(snewfun(svec2)))./Xnew(snewfun(svec2)) 0.], 'spline');
    dsC2new = griddedInterpolant( snew, [0. dsCnew(snewfun(svec2))-dsC1new(snewfun(svec2)) 0.], 'spline');
    
    % recalculate shape integral (i.e. surface area / 2pi)
    xintegral = integral(@(s) Xnew(s), 0., Lnew);

    % save arc length on initial spherical surface as function of current arc length:    
    s0new = griddedInterpolant(snew, s0(svec3));
    
    % update the tension profile according to Eqn. (D.10)
    x0 = x0fac*Lnew;    
    
    npoints = 1e+5;
    snew3 = linspace(0, Lnew, npoints); % grid for gradient calculation
    h = Lnew/(npoints-1);
    
    zetanew = griddedInterpolant(snewEu, zeta(svec3), 'spline'); % map tension profile onto updated shape
    dszetavec = gradient(zetanew(snew3), h); % recalculate tension gradient
    dszetanew = griddedInterpolant(snew3, dszetavec, 'spline');
    
end


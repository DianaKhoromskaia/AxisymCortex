function res = bc(yleft,yright, varpar,  C1, C2, C, dsC, Psi, X, Z, X0, L, zeta, dszeta, eta, etab, etac, etacb, eps1, xintegral, fext, kappa, C0, delt, xi) 
% Residual function for the bvp4c solver, specifying the boundary conditions 
% for the bvp describing the force balance of an active viscous cortex
% for a given shape at time t.
% The vector of unknowns is
% v = (v_s, d_sv_s, v_n, d_sv_n, m_s^s, t_n^s, dV(s), dX(s),
% (s0-s)_Lagrange, I1, (s0-s)_Euler)
% L -- perimeter length of half the shape

res = [ yleft(1) - 0.;    % v_s(0)=0
        yleft(4) - 0.;    % d_sv_n(0)=0
        yleft(6) - 0.;    % t_n^s(0)=0
        yleft(7) - 0.;    % dV(0)=0
        yleft(8) - 0.;    % dX(0)=0
        yleft(9) - 0.;    % (s0-s)_Lagrange(0)=0
        yleft(10) - 0.;   % I1(0)=0
        yright(1) - 0.;   % v_s(L)=0
        yright(4) - 0.;   % d_sv_n(L)=0
        yright(7) - 0.;   % t_n^s(L)=0
        yright(6) - 0.;   % dV(L)=0
        yleft(11) - 0.];  % (s0-s)_Euler(0)=0
end



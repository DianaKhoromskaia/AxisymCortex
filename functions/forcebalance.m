function [sint, Vint, P, sol, vs, dsvs, vn, dsvn, mss, tns, ds2vn, dV, dX0, Lnew, eps1new, sfun, snewfun, snewvec, sfunEu, snewfunEu, snewvecEu] = forcebalance(C1, C2, C, dsC, Psi, X, Z, X0, L, zeta, dszeta, eta, etab, etac, etacb, eps1, eps1abs, xintegral, fext, kappa, C0, P0, optode, sol, sfun, snewfun, npoints, delt, t, xi )
% This function solves the force balance for an active viscous surface at
% times t>0. The initial guess for the solver is generated from the 
% solution at the previous time step. 
% L -- perimeter length of half of the surface

%% solve force balance

solold = sol; % solution structure at the previous time step
sgrid = snewfun(solold.x); % projection of the grid at the previous time step onto current shape
sgrid = [0. sgrid(2:(end-1)) L]; % defining end and start point separately, to avoid precision issues.

parguess = [P0(1)]; % guess for parameters; P0 has two entries

solinit = bvpinit(sgrid, @yguessfun, parguess, solold, sfun); % initial guess structure for solution

sol = bvp4c( @ode, ...
    @bc, ...
    solinit,optode, C1, C2, C, dsC, Psi, X, Z, X0, L, zeta, dszeta, eta, etab, etac, etacb, eps1, xintegral, fext, kappa, C0, delt, xi);

[sint, indforward, indbackward] = unique(sol.x); % leave out any grid points which are too close to each other
Vint = sol.y(:,indforward);
derivatives = sol.yp(:,indforward);
P = sol.parameters;
dV = sol.y(7,end);
dX0 = sol.y(8,end);

%% saving new arc lengths on half-shape:
    % new arc length for Lagrange update:
Lnew = L + delt*Vint(9,end);
eps1new = eps1abs*Lnew;

sfun = griddedInterpolant( [0. sint(2:(end-1))+delt*Vint(9,2:(end-1)) Lnew], [0. sint(2:(end-1)) L], 'spline');
snewfun = griddedInterpolant( [0. sint(2:(end-1)) L], [0. sint(2:(end-1))+delt*Vint(9,2:(end-1)) Lnew], 'spline');
snewvec = [0. sint(2:(end-1))+delt*Vint(9,2:(end-1)) Lnew];

    % new arc length for Euler update: (used for active tension profile)
LnewEu = L + delt*Vint(11,end);
sfunEu = griddedInterpolant( [0. sint(2:(end-1))+delt*Vint(11,2:(end-1)) LnewEu], [0. sint(2:(end-1)) L], 'spline');
snewfunEu = griddedInterpolant( [0. sint(2:(end-1)) L], [0. sint(2:(end-1))+delt*Vint(11,2:(end-1)) LnewEu], 'spline');
snewvecEu = [0. sint(2:(end-1))+delt*Vint(11,2:(end-1)) LnewEu];

%% interpolants of flow field and derivatives, tensions, and bending moments on half-shape:
vs = griddedInterpolant(sint, Vint(1,:), 'spline');
dsvs = griddedInterpolant(sint, Vint(2,:) , 'spline');
vn = griddedInterpolant(sint, Vint(3,:) , 'spline'); 
dsvn = griddedInterpolant(sint, Vint(4,:) , 'spline');
mss = griddedInterpolant(sint, Vint(5,:), 'spline');
tns = griddedInterpolant(sint, Vint(6,:), 'spline');

ds2vn = griddedInterpolant(sint, derivatives(4,:), 'spline');

end


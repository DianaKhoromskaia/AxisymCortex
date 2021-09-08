function cortex_half_fixdt_Gaussian_cluster(factor, p, dt0)

%--------------------------------------------------------------------------
% This is the source code for simulations of the dynamics of an 
% axisymmetric viscous active cortex with an up-down symmetric profile of
% isotropic active tension. The active tension is given as a super-Gaussian. 
% The force balance at each time step is solved 
% using the ode-solver bvp4c, which yields the flow field (v_s(t), v_n(t)). 
% The equations are solved on the arc length interval s \in [0,L/2], 
% from the south pole to the equator of the shape.
% The time integration is done with the Euler method with constant time step.
%
% Input: factor -- magnitude of the active tension peak at the equator (on top of constant off-set) 
%        p      -- parameter for super-Gaussian (see below)
%        dt0    -- time step
%
% Output files: see below (lines 121-178)
%
%
% Version: 29/04/2021
% Authors: Diana Khoromskaia and Guillaume Salbreux

%--------------------------------------------------------------------------

addpath('./functions');

%% set physical paramters:

R0 = 1; % radius of initial sphere
L0 = pi*R0/2; %half-perimeter of initial sphere
L = L0;
eta = 1.; % shear 2D viscosity
etab = eta; % bulk 2D viscosity
etacb = (1e-6)*eta; % bending bulk viscosity
etac = 0.; % bending shear viscosity (not used)
kappa = 0.; % bending modulus
C0 = 0.; % spontaneous curvature
xi = 0; % effective friction with external medium

%% set numerical parameters:

%bvp-solver parameters:
epsoderel = 1e-4; % relative error tolerance
epsodeabs = 1e-6; % absolute error tolerance
optode = bvpset('RelTol', epsoderel, 'AbsTol', epsodeabs, 'stats','on','Vectorized', 'on','FJacobian',@fjac,'BCJacobian',@bcjac,'NMax', 1e+5);

% parameters for time integration:
dt = dt0; % step size
t = 0.; % start time
tmax = 5; % end time

% for adaptive time stepping: (not used)
dtmax = 1e-5; % maximal step size
tol = 1e-4; %absolute tolerance on error for time stepping
nattempt = 0; % counter for time step attempts

eps1abs = 1e-3; %cut-off at SP for time integration, as a ratio of L
eps1 = eps1abs*L; %cutoff at SP for time integratio

%% initialise spherical shape for arc length s \on [0,L0/2]
%--------------------------------------------------------------------------
% C1 = C_{\phi}^{\phi} = sin(psi)/x, circumferential curvature
% C2 = C_s^s = \partial_s psi, meridional curvature
% C = C1+C2
% X0 - centroid position of full shape
% xintegral - (surface area of full shape)/2pi
%--------------------------------------------------------------------------

z0 = 0.; %offset in z-direction
[C1, C2, C, dsC1, dsC2, dsC, Psi, X, Z, X0, xintegral0, svec1] = initialisehalfsphere(2*L0, z0);
s0 = griddedInterpolant(svec1, svec1); % arc length on initial shape
xintegral = xintegral0;
% save initial shape and arc length for plotting:
Xinit = X; 
Zinit = Z;
sinit = svec1; 

%% active tension profile 
% The profile of isotropic tension is defined as a super-Gaussian with 
% power 'pow' (pow=1 is standard Gaussian), peak height A, 
% constant off-set gamma, width parameter sigma0, and centred at x0. 

x0fac = 1; %in units of half-interval length L0/2
x0 = x0fac*L0;
A = factor;
gamma = 1; %constant background tension
pow = p/2;
sigma = 10; %sigma from Christian's formula
sigma0 = (sigma^(-1/p))/sqrt(2);

% define tension profile at t=0:
zeta = griddedInterpolant(svec1, superGaussian(svec1, x0, A, gamma, sigma0, pow));
dszeta = griddedInterpolant(svec1, dxsuperGaussian(svec1, x0, A, sigma0, pow));

% plot tension profile at t=0:
% figure(5)
% plot(svec1, zeta(svec1),'-');hold on;
% plot(svec1, dszeta(svec1),'-');
% xlabel('s');
% legend('\zeta','\partial_s \zeta');
% saveas(figure(5),strcat('zetaprofile_pbar=',num2str(pow),'_sigmabar=',num2str(sigma0),'.fig'));
% saveas(figure(5),strcat('zetaprofile_pbar=',num2str(pow),'_sigmabar=',num2str(sigma0),'.png'));

P0 = 2*gamma/R0; %Laplace pressure
P = P0;

%% define external force along along z-direction: (not used)
fextmag = 0.;
%fext = griddedInterpolant(svec1, 0.*(superGaussian(svec1, 0., -fextmag, 0., 0.05, 1)+superGaussian(svec1, fextmag, 100, 0., 0.05, 1)), 'spline');
fext = griddedInterpolant(svec1, 0.*superGaussian(svec1, 0., -fextmag, 0., 0.5, 1), 'spline');

%% sizes of different grids
npoints = 200;      %initial uniform grid for half the interval
nplot = 1000;       %for plotting and saving to file
ncompare = nplot;   %for error calculation in adaptive time stepping (not used)

%% output files:
str1 = strcat('R0=',num2str(R0),'_eta=',num2str(eta),'_zetahat=',num2str(factor),'_dt=',num2str(dt0),'_etacb=',num2str(etacb),'_zeta0=',num2str(gamma),'_Euler_pbar=',num2str(pow),'_sigmabar=',num2str(sigma0));

dir1 = strcat('results_',str1);
mkdir(dir1);
cd(dir1);

% save certain measurements on the shape at regular time steps: 
filename1 = 'observables.dat';
fileID = fopen(filename1,'w');
fprintf(fileID, ' eta=%9.8f etab=%9.8f etac=%9.8f etacb=%9.8f \n eps1abs=%e epsoderel=%e epsodeabs=%e \n dt=%e x0fac=%9.8f A=%9.8f gamma=%9.8f pow=%9.8f sigma0=%9.8f\n \n', eta, etab, etac, etacb, eps1abs, epsoderel, epsodeabs, dt, x0fac, A, gamma, pow, sigma0);
fprintf(fileID, '%s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \n', 't', 'dt', 'v_n(L)', 'P', 'r(L)', 'C2(L)','pole-pole', 'area', 'tcomp [sec]', 'L', 'dX0', 'dV','X0', 'nmesh');
formatSpec = '%6.5f \t %6.5f \t %6.5f \t %6.5f \t %6.5f \t %6.5f \t %6.5f \t %6.5f \t %6.5f \t %6.5f \t %9.8f \t %6.5f \t %6.5f \t %6.5f \n';
fprintf(fileID, formatSpec, t, dt, 0., P0, R0, 1/R0, Z(L0)-Z(0.), 2*pi*xintegral, 0., L0, 0., 0., X0, 0.);

svecuni = linspace(0., L, nplot);

% save curvatures, curvature derivatives, and tangent angle:
filename2 = 'curvatures.dat';
dlmwrite(filename2, [svecuni(1:end-1) svecuni(end)-flip(svecuni)], 'precision', '%6.5f' ,'delimiter', '\t');
dlmwrite(filename2, [C1(svecuni(1:end-1)) flip(C1(svecuni))], '-append','precision', '%6.5f' ,'delimiter', '\t');
dlmwrite(filename2, [C2(svecuni(1:end-1)) flip(C2(svecuni))], '-append','precision', '%6.5f' ,'delimiter', '\t');
dlmwrite(filename2, [dsC1(svecuni(1:end-1)) -flip(dsC1(svecuni))] , '-append','precision', '%6.5f' ,'delimiter', '\t');
dlmwrite(filename2, [dsC2(svecuni(1:end-1)) -flip(dsC2(svecuni))], '-append','precision', '%6.5f' ,'delimiter', '\t');
dlmwrite(filename2, [Psi(svecuni(1:end-1)) pi-flip(Psi(svecuni))], '-append','precision', '%6.5f' ,'delimiter', '\t');

% save x and z coordinates: 
filename3 = 'x.dat';
dlmwrite(filename3, [X(svecuni(1:end-1)) flip(X(svecuni))], 'precision', '%6.5f' ,'delimiter', '\t');
filename4 = 'z.dat';
dlmwrite(filename4, [Z(svecuni(1:end-1)) 2*Z(svecuni(end))-flip(Z(svecuni))], 'precision', '%6.5f' ,'delimiter', '\t');
filename41 = 'psi.dat';
dlmwrite(filename41, [Psi(svecuni(1:end-1)) pi-flip(Psi(svecuni))], 'precision', '%6.5f' ,'delimiter', '\t');

% save output of bvp-solver:
filename5 = 'bvpsolution.dat';
filename6 = 'vs.dat';
filename7 = 'vn.dat';
filename9 = 'zetaprofile.dat';
filename10 = 'coordinates.dat';

% save movie of the deformation:
% vidObj1 = VideoWriter('shape_sequence.mp4','MPEG-4');
% vidObj1.Quality = 100;
% vidObj1.FrameRate = 20;
% open(vidObj1);

% time step for saving into data file:
delt_data = 0.01;
n = 1.;
nsteps = 0;

% save parameters:
filename8 = 'parameters.dat';
file8 = fopen(filename8,'w');
fprintf(file8, 'R0 \t eta \t etab \t etac \t etacb \t kappa \t eps1abs \t epsoderel \t epsodeabs \t dt \t tmax \t delt_data \t A \t x0fac \t gamma \t pow \t sigma0 \t npoints \t nplot \t ncompare \t fextmag \n');
fprintf(file8, '%6.5f \t %e \t %e \t %e \t %e \t %6.5f \t %e \t %e \t %e \t %e \t %6.5f \t %e \t %6.5f \t %6.5f \t %6.5f \t %6.5f \t %6.5f \t %d \t %i \t %i \t %6.5f', R0, eta, etab, etac, etacb, kappa, eps1abs, epsoderel, epsodeabs, dt, tmax,  delt_data, A, x0fac, gamma, pow, sigma0, npoints, nplot, ncompare, fextmag);
fclose(file8);

%% simulation of viscous cortex dynamics:

while t < tmax
   
    tnew = t;
    nattempt = 0;

    %% solve force balance: 
    if t==0 % at time t=0
        tic
        [svec1, v1, P1, sol1, vs, dsvs, vn, dsvn, mss, tns, ds2vn, dV, dX0, Lnew_first, eps1_first, sfun_first, snewfun_first, snewvec, sfunEu, snewfunEu, snewvecEu] = forcebalance_t0(C1, C2, C, dsC, Psi, X, Z, X0, L, zeta, dszeta, eta, etab, etac, etacb, eps1, eps1abs, xintegral, fext, kappa, C0, P, optode, npoints, dt, xi);
        tcomp = toc; % time to compute force balance solution
    else % at times t>0
        tic
        [svec1, v1, P1, sol1, vs, dsvs, vn, dsvn, mss, tns, ds2vn, dV, dX0, Lnew_first, eps1_first, sfun_first, snewfun_first, snewvec, sfunEu, snewfunEu, snewvecEu] = forcebalance(C1, C2, C, dsC, Psi, X, Z, X0, L, zeta, dszeta, eta, etab, etac, etacb, eps1, eps1abs, xintegral, fext, kappa, C0, P, optode, sol, sfun, snewfun, npoints, dt, t, xi);
        tcomp = toc;
    end
    
    while (tnew == t)
    nattempt = nattempt+1;    
    %% evolve shape with time step dt:

    [C1_first, C2_first, C_first, X_first, Psi_first, Z_first, dsC1_first, dsC2_first, dsC_first, zeta_first, dszeta_first, xintegral_first, x0, s0new] = evolvefunctions(svec1, vs, dsvs, vn, dsvn, mss, tns, ds2vn, C1, C2, C, dsC1, dsC2, dsC, X, Psi, Z, snewfun_first, snewvec, snewfunEu, snewvecEu, t, dt, L, Lnew_first, eps1, eta, etac, etacb, gamma, x0fac, A, sigma0, pow, s0, zeta, npoints);

    scompare = linspace(0., L, ncompare);         %uniform grid on surface at time t
    sfirst = snewfun_first(scompare);               %its map on surface at time t+dt
    ssec = sfirst; %snewfun_half(scompare);                   %its map on surface at time t+dt/2 (not used)
    sthird = sfirst;%snewfun_third(snewfun_half(scompare));  %its map on surface at time t+2*(dt/2) (not used)    

    sol = sol1;
    svec = svec1;
    v = v1;
    svecuni = linspace(0., Lnew_first, nplot);     %save shape on fine uniform grid in new arc length
    svecunifb = scompare;                          %save flow field on uniform grid at time t  
    nmesh = sol.stats.nmeshpoints;                 % number of grid points in solver grid
    
    C1 = C1_first;
    C2 = C2_first;
    C = C_first;
    dsC1 = dsC1_first;
    dsC2 =  dsC2_first;
    dsC = dsC_first;
    X = X_first;
    Psi = Psi_first;
    Z = Z_first;
    s0 = s0new;
    
    X0 = X0 + dt*dX0; % the centroid is not moving in this simulation
    
    zeta = zeta_first;
    dszeta = dszeta_first; 
    xintegral = xintegral_first;
    P = P1;
    L = Lnew_first;
    eps1 = eps1_first;
    sfun = sfun_first;
    snewfun = snewfun_first;    

    tnew = t + dt;
    nsteps = nsteps + 1;
    
    end
    t = tnew;
    
    if t > n*delt_data

%     %% plot new shape:
%     
%     figure(3)
%     % plot original shape:
%     plot(Xinit(sinit), Zinit(sinit), 'b'); hold on
%     plot(-Xinit(sinit), Zinit(sinit), 'b')
%     
%     % plot shape at time t:
%     plot([X(svecuni(1:end-1)) flip(X(svecuni))], [Z(svecuni(1:end-1)) 2*Z(svecuni(end))-flip(Z(svecuni))], 'r','LineWidth',2);
%     plot(-[X(svecuni(1:end-1)) flip(X(svecuni))], [Z(svecuni(1:end-1)) 2*Z(svecuni(end))-flip(Z(svecuni))], 'r','LineWidth',2); hold off;
%     %plot(0,X0, 'ob'); 
%     axis square
%     axis([-1.5*R0 1.5*R0 X0-1.5*R0 X0+1.5*R0])
%     title(strcat('t=',num2str(t)))
        
    %% save data to file:

    dtX0 = dX0;
    dtV = dV;
    fprintf(fileID, formatSpec, t, dt, vn(L), P(1), X(L), C(L), (Z(L)-Z(0.)), 2*pi*xintegral, tcomp, L, dtX0, dtV, X0, nmesh);

    dlmwrite(filename2, [svecuni(1:end-1) 2*svecuni(end)-flip(svecuni)], '-append', 'precision', '%6.5f' ,'delimiter', '\t');
    dlmwrite(filename2, [C1(svecuni(1:end-1)) flip(C1(svecuni))], '-append', 'precision', '%6.5f' ,'delimiter', '\t');
    dlmwrite(filename2, [C2(svecuni(1:end-1)) flip(C2(svecuni))], '-append', 'precision', '%6.5f' ,'delimiter', '\t');
    dlmwrite(filename2, [dsC1(svecuni(1:end-1)) -flip(dsC1(svecuni))], '-append', 'precision', '%6.5f' ,'delimiter', '\t');
    dlmwrite(filename2, [dsC2(svecuni(1:end-1)) -flip(dsC2(svecuni))], '-append', 'precision', '%6.5f' ,'delimiter', '\t');
    dlmwrite(filename2, [Psi(svecuni(1:end-1)) pi-flip(Psi(svecuni))], '-append', 'precision', '%6.5f' ,'delimiter', '\t');

    dlmwrite(filename3, [X(svecuni(1:end-1)) flip(X(svecuni))], '-append', 'precision', '%6.5f' ,'delimiter', '\t');
    dlmwrite(filename4, [Z(svecuni(1:end-1)) 2*Z(svecuni(end))-flip(Z(svecuni))], '-append', 'precision', '%6.5f' ,'delimiter', '\t');
    dlmwrite(filename41, [Psi(svecuni(1:end-1)) pi-flip(Psi(svecuni))], '-append', 'precision', '%6.5f' ,'delimiter', '\t');

    if n==1
        dlmwrite(filename5, [svec(1:end-1) 2*svec(end)-flip(svec)], 'precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite(filename5, v, '-append', 'precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite(filename6, [vs(svecunifb(1:end-1)) -flip(vs(svecunifb))], 'precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite(filename7, [vn(svecunifb(1:end-1)) flip(vn(svecunifb))], 'precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite(filename9, [zeta(svecuni(1:end-1)) flip(zeta(svecuni))], 'precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite(filename10, scompare, 'precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite(filename10, ssec,  '-append','precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite(filename10, sfirst,  '-append','precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite(filename10, sthird,  '-append','precision', '%6.5f' ,'delimiter', '\t');
    else
        dlmwrite(filename5, [svec(1:end-1) 2*svec(end)-flip(svec)], '-append', 'precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite(filename5, v, '-append', 'precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite(filename6, [vs(svecunifb(1:end-1)) -flip(vs(svecunifb))], '-append','precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite(filename7, [vn(svecunifb(1:end-1)) flip(vn(svecunifb))], '-append','precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite(filename9, [zeta(svecuni(1:end-1)) flip(zeta(svecuni))], '-append','precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite(filename10, scompare,  '-append', 'precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite(filename10, ssec,  '-append','precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite(filename10, sfirst,  '-append','precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite(filename10, sthird,  '-append','precision', '%6.5f' ,'delimiter', '\t');
    end

    n = n + 1;
    end
end

 %close(vidObj1);
 cd ..
end

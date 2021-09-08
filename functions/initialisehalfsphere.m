function [C1, C2, C, dsC1, dsC2, dsC, Psi, X, Z, X0, xintegral, svec1] = initialisehalfsphere(L0, z0)

% This function initialises the shape descriptors as spline interpolants on
% half of the initial sphere, as functions of the arc length s \in [0,
% L0/2].

npoints = 1e+4;
svec1 = linspace(0, L0/2, npoints); %arc length variable for half of the sphere
sone = ones(1, npoints);
R0 = L0/pi;

c1vec = (pi/L0)*sone;   %C_{\phi}^{\phi} = sin(psi)/x
c2vec = c1vec;          %C_s^s = \partial_s psi
psivec = (pi/L0)*svec1;
xvec = L0*sin((pi/L0)*svec1)/pi;
zvec = R0*(sone-cos(svec1/R0));

C1 = griddedInterpolant(svec1, c1vec, 'linear');
C2 = griddedInterpolant(svec1, c2vec, 'linear');
C =  griddedInterpolant(svec1, 2*c1vec, 'linear');
dsC1 = griddedInterpolant(svec1, zeros(npoints,1),'linear');
dsC2 = griddedInterpolant(svec1, zeros(npoints,1),'linear');
dsC = griddedInterpolant(svec1, zeros(npoints,1),'linear');
Psi = griddedInterpolant(svec1, [0. psivec(2:end-1) pi/2], 'spline');
X = griddedInterpolant(svec1, [0. xvec(2:end-1) R0], 'spline');
Z = griddedInterpolant(svec1, [0. zvec(2:end-1) R0], 'spline');

xintegral = R0^2; % integral of the full shape: int_0^L x(s) ds
X0 = z0+R0; % centroid

end


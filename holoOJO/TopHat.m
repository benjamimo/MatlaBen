function U=TopHat(bbeta,r0)
% U=TopHat(X,Y,bbeta,r0)
% This function computes the phase to transform an input Gaussian beam into
% a flat-top beam.  This routine is not optimized... and indeed is SLOW!

ResR=1*1080;
% pixels
ResC=1*1920;
Ux=2*7.68E-3; % [m]
Uy=8.64E-3;

% HOLOGRAM PARAMETERS
downScale=1/1;%0.50;          % Scale factor for tests
Nx=ResC*downScale;
Ny=ResR*downScale;

% Spatial sampling at 2*fNyquist, spatial gratings
sx=linspace(-Ux/2,Ux/2,Nx);
sy=linspace(-Uy/2,Uy/2,Ny);
% [X,Y]=meshgrid(sx,sy);
% [THETA, RHO]=cart2pol(X,Y);

% Preallocate memory
PhiInt=zeros(ResR/2,ResR/2);
PhaseMap=zeros(ResR,ResR);
U=0.5*ones(ResR,ResC);

% Function to integrate
PhiFun=@(v) sqrt(1-exp(-v.^2));
counter=1;

% bbeta=20;
% r0=10.25e-4;

% Integration for different radius (really bad way to integrate :C)
for x=linspace(0,Uy,Ny/2),
    for y=linspace(0,Uy,Ny/2),
    r=sqrt(x^2+y^2);
    PhiInt(counter)=(bbeta*sqrt(pi)/4) * integral(PhiFun,0,sqrt(2)*r/r0);
    counter=counter+1;
    end
end

% Computes phase (just one quarter of the screen)
% Phase=(lambda/(2*pi*(n-1)))*mod(PhiInt,2*pi); % According to Andrews code
Phase=mod(PhiInt,2*pi); % the normalization comes after...

% To optimize we use the calculated quarter to generate the whole field
PhaseMap(1:Ny/2,1:Ny/2)=rot90(Phase,2);
PhaseMap(Ny/2+1:Ny,1:Ny/2)=fliplr(Phase);
PhaseMap(1:Ny/2,Ny/2+1:Ny)=flipud(Phase);
PhaseMap(Ny/2+1:Ny,Ny/2+1:Ny)=(Phase);

boxr=1:ResR;
boxc=(ResC-ResR)/2:(ResC-ResR)/2+ResR-1;
U(boxr,boxc) = exp(1i.*PhaseMap);
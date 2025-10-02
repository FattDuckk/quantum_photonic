%% 
% EM response of a small sphere using Mie theory - EXACT COPY
%
% wave is incident from -z -> +z
% polarized Ex

%% cleanup
clear
close all
clc

%% particle parameters
a=0.020; % sphere radius (in microns)
f=0.7; % fraction of volume occupied by metal
nc=1.5; % core index (e.g. core is silica nanoparticle)
nb=1.5; % background index (e.g. immersed in water or glass)
metal='K'; % metal of shell
[et,Lt]=feval(['eps',metal]); % check material table
Lmin = min(Lt); % ensure valid min L
Lmin = min(Lt(real(et)<0)); % only include et<0
eV = 1.24*linspace(0,1,1e3).'/Lmin; % use frequency scale (bounded)
eV = eV(eV>0); % avoid zero
L = 1.24./eV; % wavelength

%%
N=10; % number of orders: should be large enough to converge
fm=['eps',metal];
m=[nc*ones(size(L)), sqrt(feval(fm,L)), nb*ones(size(L))]; % refractive indices, inner to outer
av = a*[(1-f)^(1/3), 1]; % radii of surfaces
vol=4*pi/3*a^3;
area=pi*a^2;

% energy scale
xunit = 'eV'; xlab=xunit; fx=@(x) 1.24./x;

% section scale
ylab='C / V / k';

%% calculate spectrum
mu=ones(size(m));
C=nan([length(L),3]);
Ls=nan(size(L));
set(gca,'ColorOrderIndex',1);
tic
for k=1:length(L) % wavelength loop
    L1=L(k);
    coeffs=fmiecoeff(N, L1, av, m(k,:), mu(k,:)); % expansion coefficients
    C1=crosssection(coeffs, L1, av, m(k,:)); % calculate efficiencies
    switch ylab(1)
        case 'C'
            C1=C1*1e6; % convert to nm^2
        case 'Q'
            C1=C1/area; % "efficiency" Q, dimensionless
        case '?'
            C1=C1/vol/(2*pi/L1); % C/vol/k, dimensionless
        otherwise
            C1=C1/vol/(2*pi/L1); % C/vol/k, dimensionless
    end
    C(k,:)=C1;
    Ls(k)=fx(L(k));
    
    % plot (update)
    if (k==1)||(k==length(L))||(toc>0.2)
        plot(Ls,C)
        drawnow
        pause(0)
        tic
    end
end

% finish graph
title(sprintf('%s shell, r=%g nm, f=%g',metal, a*1e3, f))
xlabel(xlab)
ylabel(ylab)

%% find & plot peaks
[Cmax,ind]=max(C); CmaxL=Ls(ind);
hold('on')
set(gca,'ColorOrderIndex',1);
ms=12;
for k=1:3
    plot(CmaxL(k),Cmax(k),'.','MarkerSize',ms)
    text(CmaxL(k),Cmax(k),sprintf('(%g,%0.2E)',CmaxL(k),Cmax(k)))
end

% legend has to be last to avoid additional entries
legend('ext','abs','sca','Location','NorthEastOutside')
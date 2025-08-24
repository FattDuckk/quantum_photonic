%% Fresnel reflection data fitting

%% clean up
clear
close all
clc

%% starting guess goes here
n0=1+1i*1; % guess for refractive index (e.g. from Brewster)

%%
data = readtable("polarisation prac.xlsx", 'VariableNamesRange','1:1')

p_intensity = data.pintensity
p_angle = data.pangle;
p_error = data.perror;

s_intensity = data.sintensity
s_angle = data.sangle
s_error = data.serror
%% ENTER DATA HERE: column vectors, p & s same size (can fill with nan if necessary)
Ad=p_angle.'; % incidence angle (degrees)

IP=570; % p incident power
IS=570; % s incident power
RP=p_intensity.'/IP; % p reflectance
RS=s_intensity.'/IS; % s reflectance

I_err=0.029; % absolute power error (e.g. resolution of lux meter, background)
R_err=0.02; % relative power error (e.g. accuracy of lux meter, noise) 
dRP=R_err*RP+I_err/IP; % uncertainty
dRS=R_err*RS+I_err/IS;
%%

%% ENTER DATA HERE: (minimal edits; keep OG code below unchanged)

IP = 570;            % same incident beam for both polarisations
IS = 15.26;

% Use p-angle as the common grid; make everything COLUMN vectors
Ad   = p_angle(:);                 % <-- was p_angle.' (row); make it a column
RP   = p_intensity(:) / IP;
RS   = s_intensity(:) / IS;        % has NaNs where s not measured (OK)
R_err = 0.02;                      % keep your simple error model
I_err = 0.029;

% Per-point uncertainties in reflectance units (keep your additive model)
dRP  = R_err*RP + p_error(:)/IP + I_err/IP;
dRS  = R_err*RS + s_error(:)/IS + I_err/IS;

% (Optional sanity check)
% size(Ad), size(RP), size(RS), size(dRP), size(dRS)


%% plot data
hold on
set(gca,'ColorOrderIndex',1)
hold('on')
errorbar(Ad,RP,dRP,'o')
errorbar(Ad,RS,dRS,'o')  % plot the data with errorbars
xlabel('Incidence angle / degrees')
ylabel('Reflectance R')

%% fit
[x,dx]=nonlinsearch(@fresnelfunc,[real(n0);imag(n0)],[RP,RS],[dRP,dRS],Ad);
nb=x(1)+1i*x(2); % fitted refractive index
db=dx(1)+1i*dx(2); % uncertainty

%% overlay fit
A=(0:0.1:90).'; % incidence angle
Rf=fresnelfunc([real(nb),imag(nb)],A); % fitted curve
set(gca,'ColorOrderIndex',1) % reset the plot colors to match p & s data
plot(A,Rf,'-')
title([num2str(nb,'%0.2f'),' +/- ',num2str(db,'%0.2f')])
ylim([0,1])


%% overlay fit
A = (0:0.1:85).';                      % safer to stop at 85° to avoid the flick
Rf = fresnelfunc([real(nb), imag(nb)], A);  % should return N×2 [Rp, Rs]

% Plot fits colour-matched to data
set(gca,'ColorOrderIndex',1);
hPfit = plot(A, Rf(:,1), '-', 'LineWidth', 1.5);   % p fit (blue)

set(gca,'ColorOrderIndex',2);
hSfit = plot(A, Rf(:,2), '-', 'LineWidth', 1.5);   % s fit (orange)

% Redo the data plots with handles so legend works
set(gca,'ColorOrderIndex',1);
hPdata = errorbar(Ad, RP, dRP, 'o');

set(gca,'ColorOrderIndex',2);
hSdata = errorbar(Ad, RS, dRS, 'o');

% Labels and legend
xlabel('Incidence angle (degrees)');
ylabel('Reflectance R');
ylim([0,1]);
title(sprintf("Reflectance vs Incidence Angle for p- and s-Polarisation", real(nb), imag(nb), real(db), imag(db)));

legend([hPdata, hSdata, hPfit, hSfit], ...
       {'p data','s data','p fit','s fit'}, ...
       'Location','best');
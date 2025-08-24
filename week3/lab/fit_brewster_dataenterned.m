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

%% plot data
hold on
set(gca,'ColorOrderIndex',1)
hold('on')
errorbar(Ad,RP,dRP,'o')
errorbar(Ad,RS,dRS,'o')  % plot the data with errorbars
xlabel('Incidence angle / degrees')
ylabel('Reflectance R')

%% fit
%[x,dx]=nonlinsearch(@fresnelfunc,[real(n0);imag(n0)],[RP,RS],[dRP,dRS],Ad);

ydata = [RP; RS];     % 26×1 column
dy    = [dRP; dRS];   % 26×1 column
Afit  = [Ad; Ad];     % 26×1 angles (first half for p, second for s)

ff = @(p) [fresnelfunc([p(1),p(2)],Ad)(:,1);  % p reflectance
           fresnelfunc([p(1),p(2)],Ad)(:,2)]; % s reflectance

[x,dx] = nonlinsearch(ff,[real(n0);imag(n0)],ydata,dy,Afit);
nb=x(1)+1i*x(2); % fitted refractive index
db=dx(1)+1i*dx(2); % uncertainty

%% overlay fit
A=(0:0.1:90).'; % incidence angle
Rf=fresnelfunc([real(nb),imag(nb)],A); % fitted curve
set(gca,'ColorOrderIndex',1) % reset the plot colors to match p & s data
plot(A,Rf,'-')
title([num2str(nb,'%0.2f'),' +/- ',num2str(db,'%0.2f')])
ylim([0,1])



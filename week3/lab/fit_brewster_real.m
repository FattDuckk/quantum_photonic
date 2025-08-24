%% Fresnel reflection data fitting

%% clean up
clear
close all
clc

%% starting guess goes here
n0=1; % guess for refractive index (e.g. from Brewster)

%% ENTER DATA HERE: column vectors, p & s same size (can fill with nan if necessary)
Ad=[,].'; % incidence angle (degrees)

IP=; % p incident power
IS=; % s incident power
RP=[,].'/IP; % p reflectance
RS=[,].'/IS; % s reflectance

I_err=; % absolute power error (e.g. resolution of lux meter, background)
R_err=; % relative power error (e.g. accuracy of lux meter, noise) 
dRP=R_err*RP+I_err/IP; % uncertainty
dRS=R_err*RS+I_err/IS;

%% plot data
hold on
set(gca,'ColorOrderIndex',1)
hold('on')
errorbar(Ad,RP,dRP,'o')
errorbar(Ad,RS,dRS,'o') % plot the data with errorbars
xlabel('Incidence angle / degrees')
ylabel('Reflectance R')

%% fit
[x,dx]=nonlinsearch(@fresnelfunc,[real(n0)],[RP,RS],[dRP,dRS],Ad);
nb=x(1); % fitted refractive index
db=dx(1); % uncertainty

%% overlay fit
A=(0:0.1:90).'; % incidence angle
Rf=fresnelfunc([real(nb),imag(nb)],A); % fitted curve
set(gca,'ColorOrderIndex',1) % reset the plot colors to match p & s data
plot(A,Rf,'-')
title([num2str(nb,'%0.2f'),' +/- ',num2str(db,'%0.2f')])
ylim([0,1])



%% clean up
clear
close all
clc

%% inputs
metal='Cu'; % metal of sphere
Lr=[0.34, 0.7].'; % wavelength range
df=0.040; % thickness of film (microns)
nc=1.5; % input index
astep=0.005; % angle step
Ls=0.633; % selected wavelength
ottod=0.5;

Wn=1e2; % number of energies in dispersion map
Kn=5e2; % number of momenta in dispersion map

%% setup scales
tr=0:astep:90;
tc=asin(1/nc)/pi*180; % critical angle
tr=tr(tr<90);
tr=tr(tr>=tc);
%tr=tr(tr<tc*1.1);
fm=['eps',metal];
es=feval(fm,Ls);
ms=sqrt(es);
Ws=2*pi/Ls;

nreq=2; % number of modes we want to get

% k vectors
W=linspace(1./max(Lr),1./min(Lr),Wn).'*2*pi;
L=2*pi./W;
e=feval(fm,L); % permittivity of metal
m=sqrt(e); % refractive index of metal

K0=W;

K=linspace(min(K0),nc*max(K0),Kn);
K0s=Ws;

Kp=W.*sqrt(e*(1^2)./(e+1^2)); % single surface result
Kps=Ws*sqrt(es./(es+1));
Tps=asin(real(Kps)/K0s/nc)/pi*180;

K2=spp_film_solve(K0,e,df);
K2s=spp_film_solve(K0s,es,df);
T2s=asin(real(K2s)/K0s/nc)/pi*180;

tr=tr(tr<(tc+(max(T2s)-tc)*1.2));

set(gcf,'Position',get(0,'ScreenSize'))

%% contacted cover: Kretchman, SPP on far interface
d=[0, df, 0]; % physical thickness of layers
nv=[nc*ones(size(m)), m, 1*ones(size(m))]; % refractive indices - film
R=nan([length(W),length(K)]);
for w=1:length(L)
    [~,r]=ffilmstack(L(w),K(:),d,nv(w,:));
    R(w,:)=abs(r).^2; % note this isn't necessarily physical...!
end

% plot map
subplot(2,3,1)
imagesc(K,W,R)
caxis([0,1]);
axis('xy')
%zlabel('|r|^2');

h = colorbar; h.Title.String = 'R';
xlabel('k [rad \mum^{-1}]');
ylabel('\omega [rad \mum^{-1}]');
title(['Kretchman: ', num2str(df*1e3), 'nm ', metal,'; n_{in}=',num2str(nc)])

%
hold('on')
plot(real(Kp),W,'k-') % single-surface SPP
plot(nc*W,W,'c-')
plot(W,W,'c-')
%plot(2*pi/Ls*[1,nc],2*pi/Ls*[1,1],'r:') % line scan
plot(nv(w,1)*2*pi./Ls*sin(tr/180*pi),2*pi/Ls,'r.')

drawnow
pause(0)

%% Kretchman line scan
nv=[nc*ones(size(ms)), ms, 1*ones(size(ms))]; % refractive indices - film
R=nan([length(Ls),length(tr)]);
for w=1:length(Ls)
    K1=nv(w,1)*2*pi./Ls(w)*sin(tr/180*pi);
    [~,r]=ffilmstack(Ls(w),K1(:),d,nv(w,:));
    R(w,:)=abs(r).^2; % note this isn't necessarily physical...!
end

% plot map
subplot(2,3,2)
plot(tr,R,'r'); ylim([0,1]);
xlabel('\theta [^o]');
ylabel('R')

hold('on')
plot(Tps,0,'k+')

% find dips
ind1=find([0, diff(sign(diff(R))), 0]>0);
plot(tr(ind1),R(ind1),'ro')
title(sprintf('\\lambda=%d: R_{min}=%0.2e @ \\theta=%0.2f^o',Ls*1e3,R(ind1),tr(ind1)))

drawnow
pause(0)

%% fields
subplot(2,3,3)
filmplotfields(Ls,nc*K0s*sin(tr(ind1)/180*pi),d,nv);
title(sprintf('\\theta=%0.2f[^o]',tr(ind1(1))))
drawnow
pause(0)

%%-------------------------------------------------------------------------
%% gapped cover: Otto, SPP on near interface
d=[0, ottod, df, 0]; % physical thickness of layers
nv=[nc*ones(size(m)), 1*ones(size(m)), m, 1*ones(size(m))]; % refractive indices - film
R=nan([length(W),length(K)]);
for w=1:length(L)
    [~,r]=ffilmstack(L(w),K(:),d,nv(w,:));
    R(w,:)=abs(r).^2; % note this isn't necessarily physical...!
end

% plot map
subplot(2,3,4)
imagesc(K,W,R)
caxis([0,1]);
axis('xy')
%zlabel('|r|^2');

h = colorbar; h.Title.String = 'R';
xlabel('k [rad \mum^{-1}]');
ylabel('\omega [rad \mum^{-1}]');
title(['Otto: ', num2str(ottod*1e3), 'nm spacer'])

%
hold('on')
plot(real(K2(:,[1,3])),W,'k-') % film SPP
plot(nc*W,W,'c-')
plot(W,W,'c-')
%plot(2*pi/Ls*[1,nc],2*pi/Ls*[1,1],'r:') % line scan
plot(nv(w,1)*2*pi./Ls*sin(tr/180*pi),2*pi/Ls,'r.')

drawnow
pause(0)

%% Otto line scan
nv=[nc*ones(size(ms)), 1*ones(size(ms)), ms, 1*ones(size(ms))]; % refractive indices - film
R=nan([length(Ls),length(tr)]);
for w=1:length(Ls)
    K1=nv(w,1)*2*pi./Ls(w)*sin(tr/180*pi);
    [~,r]=ffilmstack(Ls(w),K1(:),d,nv(w,:));
    R(w,:)=abs(r).^2; % note this isn't necessarily physical...!
end

% plot map
subplot(2,3,5)
plot(tr,R,'r'); ylim([0,1]);
xlabel('\theta [^o]');
ylabel('R')

hold('on')
plot(T2s([1,3]),0,'k+')

% find dips
ind2=find([0, diff(sign(diff(R))), 0]>0);
plot(tr(ind2),R(ind2),'ro')

tit={};
for cp=1:length(ind2)
    tit{cp}=sprintf('R_{min}=%0.2e @ \\theta=%0.2f^o',R(ind2(cp)),tr(ind2(cp)));
end
title(tit)

drawnow
pause(0)

%% fields
subplot(4,3, 9); filmplotfields(Ls,nc*K0s*sin(tr(ind2(1))/180*pi),d,nv); title(sprintf('\\theta=%0.2f',tr(ind2(1))))
subplot(4,3,12); filmplotfields(Ls,nc*K0s*sin(tr(ind2(2))/180*pi),d,nv); title(sprintf('\\theta=%0.2f',tr(ind2(2))))

drawnow
pause(0)
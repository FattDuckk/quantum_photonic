function [kk]=spp_film_solve(k0,e,d)

% kk=spp_film_solve(k0,e,d)
%
% solve the dispersion problem for SPP film in vacuum (and single interface)
% e is the permittivity and d is the thickness

k0=k0(:); e=e(:);
k1=k0.*sqrt(e./(e+1)); % single interface
% search from k0 for the two film modes
% seems reasonably robust except at high frequency (and very low)
ki=k0;
%ki=k1;
ki=[real(ki), imag(ki)];
k2=nan(size(k0));
k3=nan(size(k0));
opt=optimset('MaxFunEvals',1e3,'MaxIter',1e3);
for kc=1:length(k0)
    [k21]=fminsearch(@(x)spp_film_func(x,k0(kc),e(kc),d, 1),ki(kc,:),opt); k2(kc,:)=k21(1)+1i*k21(2);
    [k31]=fminsearch(@(x)spp_film_func(x,k0(kc),e(kc),d,-1),ki(kc,:),opt); k3(kc,:)=k31(1)+1i*k31(2);
end

kk=[k3, k1, k2];

function m=spp_film_func(x,k0,e,d,s)

% from Economou
k=x(1)+1i*x(2);
k2=k.^2;
k02=k0.^2;
K=sqrt(k2-k02.*e);
K0=sqrt(k2-k02);
R=-(K./e)./(K0);
l=(1-R)./(1+R);
%m=abs(l).^2; % single
r=exp(-K.*d);
m=abs(l-s*r).^2;
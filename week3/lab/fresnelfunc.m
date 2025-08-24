function Rout=fresnelfunc(x,A)

% R=fresnelfunc([n,k],A)
% 
% Fresnel reflection R=[Rp,Rs]
% from air to refractive index n+ik
% incident angle A in degrees

n1=1; % incidence medium
if length(x)>1
    n2=x(1)+1i*x(2); % trasmitted medium
else
    n2=x;
end
ai=A/180*pi; % incidence angle in rad
at=asin(n1.*sin(ai)./n2); % transmitted angle
rs=(n1*cos(ai)-n2*cos(at))./(n1*cos(ai)+n2*cos(at)); % electric s
rp=(n1*cos(at)-n2*cos(ai))./(n1*cos(at)+n2*cos(ai)); % electric p
Rs=abs(rs).^2; % power s
Rp=abs(rp).^2; % power p
%% rotating analyzer output
%p=P/180*pi;
%E=(rp*sin(p)+rs*cos(p))/sqrt(2);
%% ellipsometric angles
% rho=rp./rs;
%rho=-conj(rho); % adjust for sign convention
%psi=atan(abs(rho))/pi*180;
%delta=angle(rho)/pi*180;
Rout=[Rp,Rs];


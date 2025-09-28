function b=allsbessel(N,z)

% allsbessel(N,r)
%
% generate a whole bunch of spherical bessel functions
%
% order [1:N].'
%
% returns a structure
% j & y are spherical bessels
% r means ricatti
% d means differential
% rdr is drx/r
%
% Matthew Arnold, UTS

% start with spherical bessels j & y
n=(0:N).';

% use in-built functions
s=size(z);
[nh, Z]=ndgrid(n+1/2, z(:).');
J=reshape(besselj(nh,Z),[N+1,s]);
Y=reshape(bessely(nh,Z),[N+1,s]);
f=repmat(reshape(sqrt(pi./(2*z)),[1,s]),[N+1,1]);
j=f.*J;
y=f.*Y;

% differentials
dj=dbessela(j,z);
dy=dbessela(y,z);

% kill zeroth order
j=j(2:end,:);
y=y(2:end,:);
dj=dj(2:end,:);
dy=dy(2:end,:);

% replicate radius (this may need tweaking)
s=size(z);
z=z(:);
z=reshape(z,[1, s]);
zr=repmat(z,[N, 1]);

% Ricatti-Bessels
rj=zr.*j;
ry=zr.*y;

% Ricatti-Bessel differentials
drj=j+zr.*dj;
dry=y+zr.*dy;

% rescaled Ricatti-Bessel differentials
rdrj=j./zr+dj;
rdry=y./zr+dy;

% write out arrays into a structure
b.j=j;
b.y=y;
b.h=j+1i*y;
b.dj=dj;
b.dy=dy;
b.dh=dj+1i*dy;
b.rj=rj;
b.ry=ry;
b.rh=rj+1i*ry;
b.drj=drj;
b.dry=dry;
b.drh=drj+1i*dry;
b.rdrj=rdrj;
b.rdry=rdry;
b.rdrh=rdrj+1i*rdry;
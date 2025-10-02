function [E,H,M,MU]=miefields(X,Y,Z,C,L,a,m,mu)

% [E,H,M,MU]=miefields(X,Y,Z,C,lambda,a,m,mu)
%
% calculate the mie fields
% E & H return as structures E.x, E.y, E.z etc
% M & MU show the spatial distribution of m & mu
%
% X,Y,Z are the location arrays to be calculated
% these must be the same size or scalar
% (3D arrays slow so not recommended)
%
% C contains the expansion coefficients from fmiecoeff
% 
% a is radii (inner to outer)
% m is refractive indices (inner to outer)
% mu is optional permeability (inner to outer)
%
% Matthew Arnold, UTS

% particle parameters
N=size(C,1);

if nargin<5
    mu=ones(size(m));
end

PLANE=1; % use explicit plane wave term

% derived quatities
k=m*2*pi/L;

r=sqrt(X.^2+Y.^2+Z.^2); % |r|
t=acos(Z./r); % theta
p=atan2(Y,X); % phi
s=size(r);

% order invariant functions
cp=cos(p);
sp=sin(p);
ct=cos(t);
st=sin(t);

% setup regions
A=length(a);
K=length(m);
I=cell([1, K]); % preallocate cell for indices
M=zeros(s);
MU=M;
I{1}=find(r<a(1));
M(I{1})=m(1);
for cc=2:(K-1)
    I{cc}=find((r>=a(cc-1))&(r<a(cc))); % indices for this region
    M(I{cc})=m(cc);
    MU(I{cc})=mu(cc);
end
I{end}=find(r>=a(end)); % exterior
M(I{end})=m(end);
MU(I{end})=mu(end);

% initialize fields
Er=zeros(s);
Et=zeros(s);
Ep=zeros(s);
Hr=zeros(s);
Ht=zeros(s);
Hp=zeros(s);

% rehash coefficients to make coding easier
% do in sets of four
% add zeros for y core functions
% external j=-scattered+incident
% external y=-i*scattered
if ~PLANE
    CC=[C(:,1:2), zeros([N, 2]), C(:,3:(end-2)), -C(:,(end-1):(end))+1, -i*C(:,(end-1):end)];
else
    CC=[C(:,1:2), zeros([N, 2]), C(:,3:(end-2)), -C(:,(end-1):(end)), -i*C(:,(end-1):end)]; % scattered only
end

for cc=1:K; % index
    ind=I{cc}; % select region
    R=r(ind)*k(cc); % radial factor then functions
    if ~isempty(ind) % keep freemat happy
        bes(cc)=allsbessel(N,R); % should really preallocate this
    end
end

% loop over order
n=[1:N].';
En=(i.^n).*(2*n+1)./(n.*(n+1));
% initialise Legendre function
pi1=0;

for nn=1:N
    % recursive Legendre functions
    if nn>1
        pin=((2*nn-1)*ct.*pi1-nn*pi2)/(nn-1);
    else
        pin=1;
    end
    tau=nn*ct.*pin-(nn+1)*pi1;
    pi2=pi1;
    pi1=pin;
    % angular parts
    %Mora=zeros(s);
    Mota=cp.*pin*En(nn);
    Mopa=-sp.*tau*En(nn);
    %Mera=zeros(s);
    Meta=-sp.*pin*En(nn);
    Mepa=-cp.*tau*En(nn);
    Nora=nn*(nn+1)*sp.*st.*pin*En(nn);
    Nota=sp.*tau*En(nn);
    Nopa=cp.*pin*En(nn);
    Nera=nn*(nn+1)*cp.*st.*pin*En(nn);
    Neta=cp.*tau*En(nn);
    Nepa=-sp.*pin*En(nn);
    
    % now do radial parts, looping over regions
    for cc=1:K; % index
        ind=I{cc}; % select region
        R=r(ind)*k(cc); % radial factor then functions
        
        if ~isempty(ind)
            % loop over two types of bessel waves (J & Y)
            for f=1:2
                if f==1
                    Ma=bes(cc).j(nn,:).';
                    Na=bes(cc).rdrj(nn,:).';
                else
                    Ma=bes(cc).y(nn,:).';
                    Na=bes(cc).rdry(nn,:).';
                end
                v=CC(nn,1+2*(f-1)+4*(cc-1));
                w=CC(nn,2+2*(f-1)+4*(cc-1));
                
                Nr=Ma./R;
                
                % increment fields
                Er(ind)=Er(ind)+(-i*w*Nera(ind).*Nr);
                Et(ind)=Et(ind)+(v*Mota(ind).*Ma-i*w*Neta(ind).*Na);
                Ep(ind)=Ep(ind)+(v*Mopa(ind).*Ma-i*w*Nepa(ind).*Na);
                
                Hr(ind)=Hr(ind)+(i*v*Nora(ind).*Nr)*m(cc)/mu(cc);
                Ht(ind)=Ht(ind)+(w*Meta(ind).*Ma+i*v*Nota(ind).*Na)*m(cc)/mu(cc);
                Hp(ind)=Hp(ind)+(w*Mepa(ind).*Ma+i*v*Nopa(ind).*Na)*m(cc)/mu(cc);
            end
        end
    end
    
    %[nn, N]
end

% single term plane wave
if PLANE
    pp=exp(i*k(cc)*r(ind).*ct(ind));
    Er(ind)=Er(ind)+pp.*st(ind).*cp(ind);
    Et(ind)=Et(ind)+pp.*ct(ind).*cp(ind);
    Ep(ind)=Ep(ind)-pp.*sp(ind);
    Hr(ind)=Hr(ind)+pp.*st(ind).*sp(ind)*m(cc)/mu(cc)*(-1);
    Ht(ind)=Ht(ind)+pp.*ct(ind).*sp(ind)*m(cc)/mu(cc)*(-1);
    Hp(ind)=Hp(ind)+pp.*cp(ind)*m(cc)/mu(cc)*(-1);    
end

% rescale H
mu0c=4*pi*1e-7*2.99792458e8;
Hr=-Hr/(mu0c);
Ht=-Ht/(mu0c);
Hp=-Hp/(mu0c);

% convert to rectangular
Ex=cp.*(st.*Er+ct.*Et)-sp.*Ep;
Ey=sp.*(st.*Er+ct.*Et)+cp.*Ep;
Ez=ct.*Er-st.*Et;
Hx=cp.*(st.*Hr+ct.*Ht)-sp.*Hp;
Hy=sp.*(st.*Hr+ct.*Ht)+cp.*Hp;
Hz=ct.*Hr-st.*Ht;

% pack output
E.x=Ex;
E.y=Ey;
E.z=Ez;
H.x=Hx;
H.y=Hy;
H.z=Hz;
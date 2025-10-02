function [E,H,Nx,MUx]=ffilmfields(x,ay,t,lambda,kt,d,n)

% [E,H,N,MU]=ffilmfields(x,ay,t,lambda,kt,d,n)
%
% plot fields for film stack
% E & H return as structures E.x, E.y, E.z etc
%
% x is the longitudinal coordinate as a row vector
% ay is the aperture function as a column vector
% (y is determined from kt)
% t is the transmission returned by ffilmstack
%
% for other parameters refer to ffilmstack

% set up for full space trace
a=[0, cumsum(d)]; % interfaces

% x should be row vector
x=x(:).';
if any(x<0)||any(x>a(end))
    %error('x is outside range of d')
end
X=length(x);

K=length(kt);

%
N=length(n);
sn=(real(n)>=0)-(real(n)<0); % sign of n
mu=sn*1; % permeabilities
k0=2*pi/lambda; % free space k
k=n.*k0; % wave vector magnitudes
kn=sqrt(repmat(k.^2,size(kt))-repmat(kt.^2,size(k))); % normal component
kn=kn.*repmat(sn,size(kt)); % adjust sign for negative index

% impedance constants
c=2.99792458e8;
mu0=4*pi*1e-7;
z=c*mu0*mu./n; % impedance assuming n is POSITIVE

% calculate g at all points in space
gx=zeros([K, X]);
qx=gx;
Nx=gx;
MUx=gx;
for in=1:N
    ind=find((x>=a(in))&(x<=a(in+1)));
    if ~isempty(ind) % freemat is silly
        gx(:,ind)=repmat(k(in)./kn(:,in)./z(in), [1, length(ind)]); % Hz/Ey
        qx(:,ind)=repmat(kt(:)./kn(:,in), [1, length(ind)]); % Ex/Ey
        Nx(:,ind)=n(in);
        MUx(:,ind)=mu(in);
    end
end


%% trace back in space
% from exit plane?
tx=gx;
rx=gx;
for ik=1:K
    M=[t(ik); 0]; % initialise
    M(isnan(M))=0;
    
    % working backwards
    for in=N:-1:1
        if in<N
            % loop over layers (note source->observer direction)
            p=exp(1i*kn(ik,in+1)*d(in+1)); % phase factor
            A=[1/p,0;0,p]; % phase matrix
            % field matrices F1=inv(F2)
            g1=k(in)/kn(ik,in)/z(in); % ratio H/Et for TM (p) wave
            g2=k(in+1)/kn(ik,in+1)/z(in+1); % ratio H/Et for TM (p) wave
            g=g2/g1;
            gp=1+g;
            gn=1-g;
            F=[gp, gn; gn, gp]/2; % converts + & - across interface
            
            % propagate matrix elements
            M=F*A*M;
        end
        
        ind=find((x>=a(in))&(x<=a(in+1)));
        
        if ~isempty(ind) % freemat is silly
            p=exp(1i*kn(ik,in)*(x(ind)-a(in+1))); % phase factor
            tx(ik,ind)=M(1)*p;
            rx(ik,ind)=M(2)./p;
        end
    end
end
tx((isnan(tx)))=0;
rx((isnan(rx)))=0;

% derived quantities
Ex=qx.*(-tx+rx);
Ey=tx+rx;
Hz=gx.*(tx-rx);

Ex((isnan(Ex)))=0;
Ey((isnan(Ey)))=0;
Hz((isnan(Hz)))=0;

% deal with non-normal incidence
ind=([1:K].'-floor(K/2+1));
dk=mean(kt(2:end)-kt(1:(end-1))); % freemat doesn't have diff
dy=2*pi/(dk*K);
y=ind*dy;
ind0=find(ind==0);
eib=ones(size(y));
eib=exp(1i*kt(ind0)*y);

% fourier transforms
%disp('computing FFT')
%disp('this may take a few seconds on older versions of freemat')
ak=fftshift(fft(ifftshift(ay./eib))); % aperture function in frequency space

akr=repmat(ak,size(x));

E.x=fftshift(ifft(ifftshift(Ex.*akr,1),[],1),1).*repmat(eib,[1, X]);
E.y=fftshift(ifft(ifftshift(Ey.*akr,1),[],1),1).*repmat(eib,[1, X]);
H.z=fftshift(ifft(ifftshift(Hz.*akr,1),[],1),1).*repmat(eib,[1, X]);
E.z=zeros(size(E.x));
H.y=zeros(size(H.z));
H.x=zeros(size(H.z));



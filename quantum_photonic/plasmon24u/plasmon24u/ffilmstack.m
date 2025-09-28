function [t,r,T,R]=ffilmstack(lambda,kt,d,n,tf)

% [t,r]=ffilmstack(lambda,kt,d,n,[tf])
%
% properties of a film stack
% 't' is the transmission in some sense
%
% lambda is the freespace wavelength
% kt is the lateral component of the wavevector 
%
% d is the thicknesses of the layers
% n is the complex refractive index of the layers (cover->subs)
%
% lambda or kt can be column vectors
%
% n is a row vector running from input to output
% n may also follow the dimensions of lambda
%
% tf is an optional parameter for terminations
% [0, 0] (default) 't' from one side to the other 
% [1, 0] E+ output assuming you know total tang. E input
% [1, 1] returns the eigenvalue (total->total)
% [0, 1] is not defined
%
% Matthew Arnold, UTS

k0=2*pi./lambda; % free space k

K=max(length(k0),length(kt)); % number of k vectors
N=size(n,2); % stack size

% expand quantities if necessary
k0=repmat(k0,[K,N]./size(k0));
n=repmat(n,[K,N]./size(n));
kt=repmat(kt,[K, N]./size(kt));

k=k0.*n; % wave vector magnitudes

%a=[0, cumsum(d)]; % interfaces

sn=(real(n)>=0)-(real(n)<0); % sign of n
mu=sn*1; % permeabilities (-ve n implies -ve mu)

c=2.99792458e8;
mu0=4*pi*1e-7;
z=c*mu0*mu./n; % impedance assuming n is POSITIVE

kn=sqrt(k.^2-kt.^2); % normal component
kn=kn.*sn; % change sign for negative n

g=k./kn./z; % ratio H/Et for TM (p) wave

% ok, by the time we get to here, most quantities should be [K, N]
t=zeros([K, 1]);
r=zeros([K, 1]);

% deal with input conditions
if nargin<=4
    tf=[0, 0];
end
TOTALIN=tf(1);
TOTALOUT=tf(2);
if (TOTALIN==0)&&(TOTALOUT==1)
    error('invalid input, tf=[0,1]')
end

% loops not ideal but easier than vectorising matrix operations
for ik=1:K
    M=eye(2); % initialize matrix
    
    % loop over layers (note observer->source direction)
    % i.e. backwards trace
    % generally more stable
    
    for in=N:-1:1
        p=exp(1i*kn(ik,in)*d(in)); % phase factor
        A=[1/p,0;0,p]; % phase matrix across layer
        
        % convert to total fields on output if required
        if (in==N)&&(TOTALOUT)
            M=[1, 1/g(ik,in); 1, -1/g(ik,in)]/2; % coverts total E & H to E+ & E-    
        end
        
        if in==1
            F=eye(2);
            % if required convert to total fields on input, otherwise skip
            if TOTALIN
                F=[1, 1; g(ik,in), -g(ik,in)]; % converts E+ & E- to total E & H
            end
        else
            % the general situation... transfers + & - across interface
            % this is generally better scaled than separate matrices
            gr=g(ik,in)/g(ik,in-1);
            gp=1+gr;
            gn=1-gr;
            F=[gp, gn; gn, gp]/2; % converts E+ & E- from one side to the other
        end
        
        % propagate matrix elements
        M=F*A*M;      
    end
    
    if TOTALOUT
        % find the eigenvalues
        [~,D]=eig(M);
        D=diag(D);
        t(ik)=D(1); % only need one because D(1)D(2)=1 for thin film stack? 
    else
        r(ik)=M(2,1)/M(1,1); % r
        t(ik)=1/M(1,1); % t
        % gives the + field at the output
    end
    
end

% cover any glitches
t(isnan(t))=0;

%
R=abs(r).^2; % incorrect?
T=real(g(:,end)./g(:,1)).*(abs(t).^2);

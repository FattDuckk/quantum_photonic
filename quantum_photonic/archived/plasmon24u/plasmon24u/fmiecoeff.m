function [X]=fmiecoeff(N, lambda, a, m, mu)

% X=fmiecoeff(N, lambda, a, m, [mu])
%
% solves Mie expansion coefficients
% for a multishelled spherical particle
%
% N: number of orders
% lambda: wavelength
% a: interface radii vector (from inner to outer)
% m: refractive index vector (from core to background)
% [mu]: optional permeability vector (from inner to outer)
%
% returns coefficient matrix where each row is an order
% columns are ordered:
% [core(2), shells (4 each), scattered (2)]
% in Bohren & Huffman notation for a single coat this would be
% [c d f g v w b a]
%
% Matthew Arnold, UTS, August 2007

% default values for permeability
if nargin<5
    mu=ones(size(m));
end

% derived quantities
k0=2*pi/lambda; % vacuum wave-vector
k=k0*m; % wave-vectors
n=(1:N).'; % order vector, column vector
K=length(m); % number of media
A=length(a); % number of interfaces
if K~=(A+1)
    warning('number of media and interfaces do not match')
end
L=A-1; % number of shells

% replicate
%nn=repmat(n,size(a));
a1=repmat(a.*k(1:A),size(n)); % interfaces (inner side)
a2=repmat(a.*k(2:K),size(n)); % interfaces (outer side)

% factors for H field
%h=-k./mu;
h=m./mu;
h1=h(1:A); % (inner side)
h2=h(2:K); % (outer side)

% bessel functions at interfaces

% odd theta + | odd phi - | even theta - | even phi -

b1=allsbessel(N,a1(1,:).');
b2=allsbessel(N,a2(1,:).');

M11=b1.j;
M21=b1.y;
%M31=M11+1i*M21;

M12=b2.j;
M22=b2.y;
M32=M12+1i*M22;

% odd theta + | odd phi + | even theta + | even phi -
x0=1;
%x0=k0*a(end); % scaling may help numerical results

N11=b1.rdrj*x0;
N21=b1.rdry*x0;
%N31=N11+1i*N21;

N12=b2.rdrj*x0;
N22=b2.rdry*x0;
N32=N12+1i*N22;

% field matching
i2=(1:L);
i1=(2:A);

C=4+4*L; % number of coefficients
X=nan([N, C]); % pre-allocate output
% [c d f1 g1 v1 w1 b a]

% assume no coupling between orders
% loop over order

for q=n.'
    B=zeros([C, 1]); % pre-allocate RHS
    M=zeros([C, C]); % pre-allocate matrix
    
    % E theta (real)
    M(0*A+1,1)=M11(q,1); % [c] (inner) core
    if L>0
        M(sub2ind(size(M),0*A+i1,-1+4*i2))=M11(q,i1); % [f] inner shells
        M(sub2ind(size(M),0*A+i1, 1+4*i2))=M21(q,i1); % [v] inner shells
        M(sub2ind(size(M),0*A+i2,-1+4*i2))=-M12(q,i2); % [f] outer shells
        M(sub2ind(size(M),0*A+i2, 1+4*i2))=-M22(q,i2); % [v] outer shells
    end
    M(0*A+A,C-2+1)=--M32(q,end); % [b] (outer) scattered
    B(0*A+A)=M12(q,end); % (outer) incident
    
    % E theta (imag)
    M(1*A+1,2)=-N11(q,1); % [d] (inner) core
    if L>0
        M(sub2ind(size(M),1*A+i1, 0+4*i2))=-N11(q,i1); % [g] inner shells
        M(sub2ind(size(M),1*A+i1, 2+4*i2))=-N21(q,i1); % [w] inner shells
        M(sub2ind(size(M),1*A+i2, 0+4*i2))=--N12(q,i2); % [g] outer shells
        M(sub2ind(size(M),1*A+i2, 2+4*i2))=--N22(q,i2); % [w] outer shells
    end
    M(1*A+A,C-2+2)=-N32(q,end); % [a] (outer) scattered
    B(1*A+A)=-N12(q,end); % (outer) incident
    
    % H theta (imag)
    M(2*A+1,1)=h1(1)*N11(q,1); % [c] (inner) core
    if L>0
        M(sub2ind(size(M),2*A+i1,-1+4*i2))=h1(i1).*N11(q,i1); % [f] inner shells
        M(sub2ind(size(M),2*A+i1, 1+4*i2))=h1(i1).*N21(q,i1); % [v] inner shells
        M(sub2ind(size(M),2*A+i2,-1+4*i2))=-h2(i2).*N12(q,i2); % [f] outer shells
        M(sub2ind(size(M),2*A+i2, 1+4*i2))=-h2(i2).*N22(q,i2); % [v] outer shells
    end
    M(2*A+A,C-2+1)=--h2(end).*N32(q,end); % [b] (outer) scattered
    B(2*A+A)=h2(end)*N12(q,end); % (outer) incident
    
    % H theta (real)
    M(3*A+1,2)=-h1(1).*M11(q,1); % [d] (inner) core
    if L>0
        M(sub2ind(size(M),3*A+i1, 0+4*i2))=-h1(i1).*M11(q,i1); % [g] inner shells
        M(sub2ind(size(M),3*A+i1, 2+4*i2))=-h1(i1).*M21(q,i1); % [w] inner shells
        M(sub2ind(size(M),3*A+i2, 0+4*i2))=--h2(i2).*M12(q,i2); % [g] outer shells
        M(sub2ind(size(M),3*A+i2, 2+4*i2))=--h2(i2).*M22(q,i2); % [w] outer shells
    end
    M(3*A+A,C-2+2)=-h2(end).*M32(q,end); % [a] (outer) scattered
    B(3*A+A)=-h2(end).*M12(q,end); % (outer) incident
    
    rX=rank(M);
    if (q>rX)
        %[q, rX, log10(cond(M))]
        break
    end
    %
    if rcond(M)<eps
        break
    end
    %[q, rX]
    X(q,:)=(M\B).'; % solve for coefficients
    %X(q,:)=(inv(M)*B).';
    
end

X=X(all(~isnan(X),2),:);




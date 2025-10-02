function [S1,S2,S11,S12,S33,S44]=mieangles(C, theta)

% [S1,S2,S11,S12,S33,S44]=mieangles(C, theta)
%
% calculate angular scattering parameters
%
% assumes C from fmiecoeff
%
% ref: Bohren & Huffman
%
% Matthew Arnold, UTS
%

N=size(C,1);
an=C(:,end);
bn=C(:,end-1);

S1=0;
S2=0;

% order invariant functions
ct=cos(theta);
st=sin(theta);

% loop over order
n=[1:N].';
En=(2*n+1)./(n.*(n+1));
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
    
    S1=S1+En(nn)*(an(nn)*pin+bn(nn)*tau);
    S2=S2+En(nn)*(an(nn)*tau+bn(nn)*pin);
    
end

S11=(abs(S1).^2+abs(S2).^2)/2;
S12=(abs(S2).^2-abs(S1).^2)/2;
S33=(conj(S2).*S1+S2.*conj(S1))/2;
S44=(conj(S2).*S1-S2.*conj(S1))/2;

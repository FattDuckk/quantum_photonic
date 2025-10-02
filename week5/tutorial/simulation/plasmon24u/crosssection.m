function C=crosssection(X, w, a, m)

% C=crosssection(X, lambda, a, m)
%
% calculates scattering cross-sections for a spherical particle
%
% [C_ext, C_abs, C_sca]
%
% assuming we have expansion coefficients in X
% where the row is order
% columns go [... b a]
%
% assume we have k=2*pi*m/lambda
%
% assume we have size of particle 'a'
%
% Matthew Arnold, UTS

k=2*pi*m(end)/w;

%x=k(end)*a(end); % size parameter

n=(1:size(X,1)).'; % order vector
% extract scattering coefficients
an=X(:,end);
bn=X(:,end-1);

kb=k(end); % k in background
C_sca=2*pi/(kb.^2)*sum((2*n+1).*(abs(an).^2+abs(bn).^2));
C_ext=2*pi/(kb.^2)*sum((2*n+1).*real(an+bn));
C_abs=C_ext-C_sca;
%G=pi*(a(end)^2); % geometric

C_sca=C_sca;
C_ext=C_ext;
C_abs=C_abs;
%Q_b=(abs(sum((2*n+1).*((-1).^n).*(an-bn)))/x)^2;

%Q_sca_small=8/3*(x^4)*(abs((m^2-1)/(m^2+2))^2)

C=real([C_ext, C_abs, C_sca]);
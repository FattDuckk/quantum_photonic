function [p,up,Dn,y0]=nonlinreg(funfcn,p,ya,dya,varargin)

% [p,up,Dn,y0]=nonlinreg(funfcn,p,ya,dya,varargin)
%
% non-linear regression
%
% funfcn is model function
% p are parameters, with uncertainty up
% ya is data, with uncertainty dya
% all input data is contained in varargin
%
% Dn is fit norm
% y0 is best fit

ind=(~isnan(ya))&(dya>eps); % will remove all invalid data
ya=ya(ind);
dya=dya(ind);

N=length(ya);
M=length(p);
tol=max(0.5*sqrt(N-M),1e-6); % target tolerance for fit
itmax=1e3; % max iterations (e.g. 1e3)

dp=sqrt(eps); % displacement for derivative estimates
dp=1e-14;
Dn=inf; % initial norm
W=diag(1./dya); % weighting

% iterate
it=0;
while (Dn>tol)
    % calculate the Jacobian
    y0=feval(funfcn,p,varargin{:}); % fitted function @ current point
    y0=y0(ind);
    J=zeros([N,M]); % place-holder
    for ip=1:M % for each parameter calc derivative
        pp=p;
        pp(ip)=pp(ip)+dp; % perturbed single parameter
        yp=feval(funfcn,pp,varargin{:}); % perturbed output
        yp=yp(:);
        J(:,ip)=(yp(ind)-y0)/dp;
    end
    Q=W*J; % weighted
    iC=Q.'*Q;
    C=inv(iC); % auxiliary matrix
    % check for instability
    if cond(iC)>1e7
        warning('fit failed - try bigger errors, different starting values, less parameters')
        break
    end
    Dy=ya-y0; % residual
    Dw=W*Dy; % weighted
    Dn=norm(Dw); % norm
    p=C*(Q.'*Dw)+p; % update params, steepest descent
    it=it+1;
    if it>=itmax
        warning('reached max iteration - fit incomplete');
        break
    end
end
up=sqrt(diag(C)); % uncertainty in params
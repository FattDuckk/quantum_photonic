function [p,up,Dn,y0]=nonlinsearch(funfcn,p,ya,dya,varargin)

% [p,up,Dn,y0]=nonlinsearch(funfcn,p,ya,dya,varargin)
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
tol=max(sqrt(N-M)/3,1e-6); % target tolerance for fit
% target chi-squared value
% divisor is to account for unc = 3 std

W=diag(1./dya); % weighting

srs=struct('type','()','subs',{{ind}});

ff=@(p1) feval(funfcn,p1,varargin{:}); % the function to fit to
Df=@(p1) norm(W*(ya-subsref(ff(p1),srs))); % norm with masking

%% fit
[p, Dn]=fminsearch(Df, p); % use fminsearch to do the non-linear search
Dn/tol
if Dn>tol
    warning('inconsistent fit - try increasing uncertainties or remove bad data')
end
y0=ff(p);

%% find uncertainties
% use non-linear search to find (diag) uncertainties
up=nan(size(p));
for ip=1:M
    % vary each parameter in turn
    % find p for which Dn*(Df-Dn)*2 ~ 1
    sr1=struct('type','()','subs',{{ip}});
    pm1=@(p1) subsasgn(p,sr1,p1); % p with single value modified
    Df1=@(p1) Dn*(Df(pm1(p1))-Dn)*2-1; % uncertainty function -> 0
    pp1=fzero(Df1,p(ip));
    dp1=pp1-p(ip);
    % trying to find opposite side seems difficult - leave as is
    up(ip)=abs(dp1);
end


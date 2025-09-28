function filmplotfields(L,K,d,nv)

%% field plotting section
%%% attempt to plot the fields
w=2*2*pi/real(K); % lateral size of plot
K0=2*pi/L;
Kn=sqrt(K0^2-K^2);

X=500; % number of long. samples
%x=linspace(0,sum(d),X);
df=d(end-1);
de=2/imag(Kn);
d(end)=d(end)+de; % extend the endpoint for plotting
cd=cumsum(d);
x=linspace(cd(end-2)-de,cd(end),X);
%x=linspace(0,cd(end),X);

% setup k vectors
nK=200;
ind=((1:nK).'-floor(nK/2+1));
dy=w/nK;
y=ind*dy;
kt=2*pi*ind/(dy*nK); % transverse wave vector (i.e spatial frequency)
kt=kt+K;

% incident field
ay=exp(1i*K*y);
% refractive indices defined above
%% calculate
t=ffilmstack(L,kt,d,nv); % calculate transmission
[E,H]=ffilmfields(x,ay,t,L,kt,d,nv); % calculate fields

%%
S=poynting(E,H);
e=fieldnorm(E);
h=fieldnorm(H);
s=fieldnorm(S);

a=[0, cumsum(d)];

%%
imagesc(x,y,real(E.x)) % x is normal so looks like charge
%imagesc(x,y,real(E.y))
%rgb=hsv2rgb(cat(3,(angle(E.x)/pi+1)/2,ones(size(e)),abs(E.x)/max(abs(E.x(:))))); % complex E.x
%rgb=hsv2rgb(cat(3,(atan2(real(E.y),real(E.x))/pi+1)/2,ones(size(e)),e/max(e(:)))); % field direction
%image(x,y,rgb)
axis('tight')
aa=axis;
colormap(gca,cmapwave(256))
%colormap(gca,'hsv')
h = colorbar; h.Title.String = 'E_x';

%colormap(jet)
hold on
plot([a; a], repmat([min(y); max(y)],[1, length(a)]),'k-') % plot the boundaries
axis(aa); % reset boundaries
%title('Ex');

%% streamline power
% boring
if 1
    h=streamslice(x,y,S.x,S.y,0.25); % power
    for k=1:length(h)
        set(h,'Color','k');
    end
end

%% streamline field
if 1
    % messy and a bit uneven
    h=streamslice(x,y,real(E.x),real(E.y),0.25,'noarrows'); 
    for k=1:length(h)
        set(h,'Color','w');
        set(h,'LineStyle','--');
    end
end

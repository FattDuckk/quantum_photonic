function y=wave(n)
   
   % y=wave(n)
   %
   % colormap
   
   if nargin<1
      n=64; % number of colors
   end
   h=floor(n/2);
   a=[h:-1:1].'/h;
   b=[1:h].'/h;
   if (n-2*h)>0
      b=[0;b];
   end
   ua=ones(size(a));
   ub=ones(size(b));
   %y=[a*[1,0,0];b*[0,0,1]];
   y=[[(1-a); ub],[(1-a);(1-b)],[ua;(1-b)]];
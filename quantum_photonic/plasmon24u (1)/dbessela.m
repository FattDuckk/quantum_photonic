function dy=dbessela(y,x)

% dbessel(y,x)
%
% differentiate spherical bessel by recursion
%
% assume that n=[0:(N-1)].'
%
% Matthew Arnold, UTS

s=size(x);
x=x(:).';
N=size(y,1);
y=reshape(y,[N, length(x)]);
dy=zeros(size(y));

v=[0:(N-1)].';

% not the most effective but it works
for k=2:N;
    n=v(k);
    dy(k,:)=y(k-1,:)-(n+1)./x.*y(k,:);
end

%dy=dy(2:(N-1),:);

dy=reshape(dy,[size(dy,1), s]);
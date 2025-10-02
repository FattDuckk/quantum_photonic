function S=poynting(E,H)
   
   % S=poynting(E,H)
   %
   % calculates time averaged power flow
   %
   % all quantities assumed to be structures S.x S.y S.z etc
   %
   % Matthew Arnold, UTS
   
S.x=real(E.y.*conj(H.z)-E.z.*conj(H.y))/2;
S.y=real(E.z.*conj(H.x)-E.x.*conj(H.z))/2;
S.z=real(E.x.*conj(H.y)-E.y.*conj(H.x))/2;
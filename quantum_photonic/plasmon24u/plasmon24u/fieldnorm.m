function f=fieldnorm(F)
   
   % function f=fieldnorm(F)
   % 
   % calculates magnitude f of vector field F
   % assuming structure F.x F.y F.z etc
   %
   % Matthew Arnold, UTS
   
   f=sqrt(abs(F.x).^2+abs(F.y).^2+abs(F.z).^2);
   
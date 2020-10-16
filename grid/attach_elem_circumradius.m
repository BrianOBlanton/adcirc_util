function fgsout=attach_elem_circumradius(fgsin)
%ATTACH_ELEM_CIRCUMRADIUS attach element circumradius to a 
% fem_grid_struct
% Call as: fgsout=attach_elem_circumradius(fgsin);

if nargin==0 && nargout==0
   disp('fgsout=attach_elem_circumradius(fgsin);')
   return
end

% Copy incoming fem_grid_struct for operations
fgsout=fgsin;

S=abs(fgsin.A_cart+1i*fgsin.B_cart);
a=S(:,1);
b=S(:,2);
c=S(:,3);
numer=(a.*b.*c);
denom=sqrt((a + b + c).*(b + c - a).*(c + a - b).*(a + b - c));
fgsout.CR = numer./denom;


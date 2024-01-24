function U = rotphi(ang,phi,spin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%U = rotphi(ang,phi,spin);
%ang - Rotation angle
%phi -define the rotation axis on XY plane (cos(phi),sin(phi),0) 
%spin- spin value (1/2,1,3/2,2, ... ), if omitted default is 1/2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin == 2; spin = 1/2; end

[Ix, Iy, Iz] = mat_ixyz(spin);

U = expm(-1i * (cos(phi) * Ix + sin(phi) * Iy ) * ang);


end
 

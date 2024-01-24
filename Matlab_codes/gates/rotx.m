function U = rotx(ang,spin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%U = rotx(ang,phi,spin);
%ang - Rotation angle (around x direction)
%spin- spin value (1/2,1,3/2,2, ... ), if omitted default is 1/2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin == 1; spin = 1/2; end

[Ix, Iy, Iz] = mat_ixyz(spin);

U = expm(-i * Ix * ang);

end
 

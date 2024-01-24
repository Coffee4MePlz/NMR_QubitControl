function U = rotz(ang,spin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%U = rotz(ang,phi,spin);
%ang - Rotation angle (around z direction)
%spin- spin value (1/2,1,3/2,2, ... ), if omitted default is 1/2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin == 1; spin = 1/2; end

[Ix, Iy, Iz] = mat_ixyz(spin);

U = expm(-i * Iz * ang);

end
 

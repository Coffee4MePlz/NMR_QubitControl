function [Ix, Iy, Iz] = mat_ixyz(In);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Alexandre M. de Souza
% [Ix, Iy, Iz] = Mat_Ixyz(In);  
%
%  In - Spin value (Ex: 1/2, 3/2, 1, ...)
%  [Ix, Iy, Iz] spin operators (units of  h/ 2*pi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


mI = In:-1:-In;  

%---------------------- Calculating the marix Ix, Iy and Iz ------------------------------

for lin=1:1:length(mI)
   for col=1:1:length(mI)
      
      if lin == col dmm = 1; else dmm = 0; end;
      
      if lin == col-1 dmp = 1; else dmp = 0; end;
      
      if lin == col+1 dpm = 1; else dpm = 0; end;
      
  		Ip = sqrt((In - mI(col)) * (In + mI(col) + 1)) * dmp;
      
      Im = sqrt((In + mI(col)) * (In - mI(col) + 1)) * dpm;
          
      
      Ix(lin, col) = (Ip + Im)/2; 
      
      Iy(lin, col) = (Ip - Im)/(2*i); 
      
      Iz(lin, col) = mI(col) * dmm;  

   end;
end;

%-----------------------------------------------------------------------------------------

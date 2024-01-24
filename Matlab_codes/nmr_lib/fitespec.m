function spec = fitespec(spec,freq,T2,range,frrang);
%Fit N peaks
%
%function -> spec = fitespec(spec,freq,T2,range,fignum,frrang);
%spec ->spectrum structure 
%freq  ->frequencies where the picks are located
%frangr ->range of frequency allowed in the fit
%T2 -> initial guess for T2
%range -> range for freq and T2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create frequency axis
 dt = 1/spec.sw;
 SI = length(spec.esp);
 fr = (1/dt)*(0:SI-1)/SI - spec.sw/2; 

 x1 = round( (frrang(1) +  spec.sw/2 )*dt*SI) +1;
 x2 = round( (frrang(2) +  spec.sw/2 )*dt*SI) +1;

 opt=optimset('TolFun',1e-20,'TolX',1e-16,'MaxIter',15000);
 
 fr = fr(x1:x2);
 esp = spec.esp(x1:x2);
optimset('TolFun',1e-20,'MaxFunEvals',5e10);
 %Funlist is a function with 64 lorentzians
for k=1:length(freq)
    funlist{2*k-1} = @(c,x) (1/(2*pi*c(1)))./((x - (freq(k) + c(2))).^2 + (1/(2*pi*c(1)))^2);
    funlist{2*k}   = @(c,x) (x- (freq(k) + c(2)))./((x - (freq(k) + c(2))).^2 + (1/(2*pi*c(1)))^2);
end
 funlist{2*k+1}   = @(c,x) 1;
[spec.INLP spec.ILP] = fminspleas(funlist,[T2 0],fr,real(esp),...
	[T2-range(1) -range(2)],[T2+range(1) range(2)],[],opt);
spec.fr = freq;
function  struc = decompulse(shape,name,pw,targ,ang,typs);
% usage --> pulso{1} = decompulse({{'isech200.RF'}},'F1180',0.3e-3,1,180,{{'n'}});

% shps := {{'pulse type.RF'}} cell array % in the format of .RF [phase,amplitude, timestep]
% name = name of the output structure, ex: 'gaussian', 'F1180', etc
% pw := float % pulse width,     %  targ := vector of int % target qubits
% ang := float / maybe an array (?) % angle of rotation
% typs := cell array {{'n'}} for normal (?) or {{'g'}} for gradient (?)
% nchan = number of channels

% optimization settings (?) (seems to not be used anywhere) and are somewhat irrealistic...
opt = optimset('TolX',1e-20,'TolFun',1e-20,'MaxFunEvals',5e10);

global mol spectro

% storing variables in a structure for output
struc.name = name;
npul = length(shape);
struc.forma = shape;
struc.targ = targ;
struc.pw = pw;
struc.ang = ang;


nchan = 0;

for m=1:npul
     if any(targ(m) == spectro.chanel{1});  nchan =nchan+1;  spintarg{nchan} = spectro.chanel{1}; struc.chan{nchan} = [1]; end 
     if any(targ(m) == spectro.chanel{2});  nchan =nchan+1;  spintarg{nchan} = spectro.chanel{2}; struc.chan{nchan} = [2]; end 
end

Uid = 1;
n=1;

for k=1:nchan
        flag(k) = 0;
        clear shp targ2 ang2 phase2;
        
    for m =1:length(shape{k})
		    shp{m} = load([shape{k}{m}]);
            if typs{k}{m} == 'n';
            flag(k) = 1;
            ang2(m) = ang(n);
            phase2(m) = 0;
            targ2(m) = targ(n);
            Uid = multigate(mol.nspin,targ2(m),rotx(ang2(m)*pi/180))*Uid; 
                %U = multigate(5,[2 4],Uop) -> U is a 5 qubit unitary where Uop is applied to qubits 2 and 4
            typs2{m} = 'n';
            end
            n =n+1;
        
        
        if flag(k)
        [wr{k} ph{k} aa np dt] = composepul(shp,pw,ang2,phase2,targ2,typs2);
         struc.power(k) = aa;
	    else
		    struc.power(k) = max(shp{m}(:,2))/(2*pi);
        end
           
        struc.typs = typs;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
struc.za = zeros(1,mol.nspin);struc.zp  = struc.za;
struc.zza = zeros(mol.nspin,mol.nspin);struc.zzp  = struc.zza;

flag =flag(flag==1);
if any(flag)
    for m=1:length(flag) % length flag??? how can this be different than 1 ?
    % defines Ix{m} and Iy{m} as the sum of all Ix's and Iy's
        Ix{m} = 0; Iy{m} = 0; % these are empty cells...
    
         for n=1:length(spintarg{m})
                Ix{m} = Ix{m} + mol.Ix{spintarg{m}(n)};
                Iy{m} = Iy{m} + mol.Iy{spintarg{m}(n)};
        end
         
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Ut = 1; Urf =1; Uz =1;
    
    
    for k=1:np 
             Hrf{k} = 0;
             
             for m=1: nchan
                 
             Hrf{k} =   Hrf{k} +  wr{m}(k)*(cos(ph{m}(k))*Ix{m} + sin(ph{m}(k))*Iy{m});
    
             end
             
             Ut = expm(-1i*(Hrf{k} +  mol.Hint + mol.Hzee)*dt)*Ut;
             
             Uz = expm(-1i*(Hrf{k} +  mol.Hzee)*dt)*Uz;
             
             Urf = expm(-1i*(Hrf{k})*dt)*Urf;
             
    end
    Ut = expm(-1i*(mol.Hint + mol.Hzee)*spectro.rof2*1e-6)*Ut*...
         expm(-1i*(mol.Hint + mol.Hzee)*spectro.rof1*1e-6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Z
        Ur1 = expm(-i*(mol.Hzee)*spectro.rof1*1e-6); 
        Ur2 = expm(-i*(mol.Hzee)*spectro.rof2*1e-6);
        
        %comparing the difference between Uz (zeeman + rf) delayed and Uid (multigate unitary) generated above
        fidelz = @(theta)  1 - fidelmat(Ur2*Uz*Ur1, errz(theta(1:mol.nspin),mol) * ...
            Uid *errz(theta(mol.nspin+1:2*mol.nspin),mol));
    
        x = fminsearch(fidelz,[pi*pw*mol.dq pi*pw*mol.dq]);
    
        struc.zp =mod(x(1:mol.nspin)*180/pi,360);
        struc.za =mod(x(mol.nspin+1:2*mol.nspin)*180/pi,360);
        Urza = errz(x(mol.nspin+1:2*mol.nspin),mol);
        Urzp = errz(x(1:mol.nspin),mol);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ZZ
    
    Urzza =1;
    Urzzp =1;
    
    struc.zzp(1:mol.nspin,1:mol.nspin) = 0;
    struc.zza(1:mol.nspin,1:mol.nspin) = 0;
    
    for k=1:mol.nspin
        
        for n=k+1:mol.nspin
        
        Uz = 1;
        zzaux = 2*pi*mol.j(k,n)*mol.Iz{k}*mol.Iz{n};
        Uz = expm(-1i*(zzaux)*spectro.rof1*1e-6);
        for m=1:np  Uz = expm(-1i*(Hrf{m} +  zzaux)*dt)*Uz; end
        Uz = expm(-1i*(zzaux)*spectro.rof2*1e-6)*Uz;
        
        fidelz = @(theta) 1 - ...
        fidelmat(Uz, zzgate(mol.nspin,theta(1),k,n) * Urf * ...
        zzgate(mol.nspin,theta(2),k,n) );
    
        x = fminsearch(fidelz,[0 0]);
    
        struc.zzp(k,n) = x(1)*180/pi;
        
        struc.zza(k,n) = x(2)*180/pi;   
       
        Urzzp =  zzgate(mol.nspin,x(1),k,n)*Urzzp; 
        Urzza =  zzgate(mol.nspin,x(2),k,n)*Urzza; 
        
        
        end
       
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    struc.f = fidelmat(Ut,Urzzp*Urzp*Uid*Urza*Urzza);
    
    disp(struc.f)

end
end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   function S = errz(theta,mol); 
   
   S = 1;
   
   for k=1:mol.nspin
         S =  multigate(mol.nspin,k,rotz(theta(k)))*S; 
   end
   end
    




    



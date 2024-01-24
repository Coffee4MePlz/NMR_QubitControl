
function  strucout = compilaSeqB(seq,pulse,seq2);
%strucout = compilaSeqB(seq,pulse);
%Compile a sequence "seq" 
%
clear  aux struc;
global mol spectro aux struc;

aux.opt = optimset('TolFun',1e-14,'MaxFunEvals',5e12,'TolX',1e-14);
aux.N = mol.nspin;          
aux.count = 0;
aux.count2 = 1;
aux.ttot = 0; 
aux.ttotp = 0;
aux.flag = 0;
aux.listspin(1:aux.N) = 0;  
aux.listvirtua(1:aux.N) = 0;  
aux.phasevirtua(1:aux.N) = 0;  
aux.Uid{1} = eye(2^aux.N);     
aux.Usim{1} = eye(2^aux.N);      
aux.errozz(1:aux.N,1:aux.N) = 0; 
aux.p = pulse;               
aux.atarg = [];
aux.atargr = []; 
aux.refoc = [];
aux.erroz(1:aux.N) = 0;
aux.list0=aux.listspin;
eval(seq); 

finseq;

disper('final')

strucout=writeoutput(struc);
strucout = writepul(strucout,seq,seq2);
for k=1:(aux.count2); 
disp(['fidelity = ' mat2str(fidelmat(aux.Usim{k},aux.Uid{k}))])
end

disp(['Total time = ' mat2str(aux.ttot*1e3) ' ms'])
disp(['Total time with pulses = ' mat2str(aux.ttotp*1e3) ' ms'])
clear global struc aux


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function appulse(nome,phase,optin);
global aux struc;

    pul=findpul(nome,aux.p);   
    targ = aux.p{pul}.targ;   ang = aux.p{pul}.ang;  pw = aux.p{pul}.pw;
    zzpos = aux.p{pul}.zzp;   zzpre = aux.p{pul}.zza;
    zpre = aux.p{pul}.za;     zpos = aux.p{pul}.zp;
    ntarg = length(targ);
    
%--------------------------------------------------------------------------
 
if aux.flag == 0;
    aux.errozz = reseta(aux.errozz + propagv(zzpre));
    aux.erroz= mod( reseta(aux.erroz + propagv(zpre)),360);
else
   optdelay(zzpre,zpre);
   aux.flag = 0;
end

%--------------------------------------------------------------------------


aux.count = aux.count +1; 
struc{aux.count}.tip = 'rotzp';
struc{aux.count}.desc = 0*(1:aux.N);
struc{aux.count}.desc(targ) = -aux.erroz(targ);
struc{aux.count}.desc  = propagv(struc{aux.count}.desc);
aux.erroz(targ) = 0;

for m=1:ntarg
   for k=1:aux.N; aux.Usim{aux.count2} = zzgate(aux.N,(pi/180)*(aux.errozz(targ(m),k) +...
        aux.errozz(k,targ(m))),targ(m),k)*aux.Usim{aux.count2}; 
   end
        aux.errozz(targ(m),:) = 0; aux.errozz(:,targ(m)) = 0; aux.listspin(targ(m)) = 0;
end

if nargin == 3; aux.listspin(optin) = 1; end

%--------------------------------------------------------------------------
for m=1:ntarg

    phasep(m) = phase(m); rz(m) = 0;

    if aux.listvirtua(targ(m)) == 1; 
        rz(m) = mod(2*(phase(m) - aux.phasevirtua(targ(m))),360);
        phasep(m) = 2*aux.phasevirtua(targ(m)) - phasep(m) - 180;
        aux.listvirtua(targ(m)) = 0;
        aux.phasevirtua(targ(m)) = 0;        
    end
    phasep(m) = mod(phasep(m),360);
end

%--------------------------------------------------------------------------

for m=1:ntarg
     V = multigate(aux.N,targ(m),rotphi(aux.p{pul}.ang(m)*pi/180,phase(m)*pi/180));
     aux.Usim{aux.count2} = V*aux.Usim{aux.count2};
     aux.Uid{aux.count2} = V*aux.Uid{aux.count2};
end


aux.count = aux.count +1; 
struc{aux.count}.tip = 'pulse'; 
struc{aux.count}.desc = aux.p{pul}; 
struc{aux.count}.phase = phasep; 
struc{aux.count}.list = aux.listspin; 

aux.count = aux.count +1; 
struc{aux.count}.tip = 'rotzp';
struc{aux.count}.desc = 0*(1:aux.N);
struc{aux.count}.desc(targ) = rz;

%--------------------------------------------------------------------------

aux.errozz = mod(reseta(aux.errozz + propagv(zzpos)),720);
aux.erroz= mod( reseta(aux.erroz + propagv(zpos)),360);

aux.ttot = aux.ttot + pw;  
aux.ttotp = aux.ttotp  + pw;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function optdelay(zzang,zang)
global mol spectro aux struc;

zzpre{1} = 0; zzpos{1} = 0;  zpre{1} = 0;  zpos{1} = 0;  ntarg{1} = 0; pw{1} = 0;

for k=1:length(aux.refoc);
     pul=findpul(aux.refoc{k}.nome,aux.p);   
     targ{k} = aux.p{pul}.targ;  ang{k} = aux.p{pul}.ang;  pw{k} = aux.p{pul}.pw;     
     zzpos{k} = aux.p{pul}.zzp;  zzpre{k} = aux.p{pul}.zza;
     zpre{k} = aux.p{pul}.za;    zpos{k} = aux.p{pul}.zp;  
     phase{k} = aux.refoc{k}.phase; ntarg{k}= length(targ{k}); 
end

%--------------------------------------------------------------------------
     par{1}.zz = mod(reseta(aux.errozz + propagv(zzpre{1})),720);
     par{1}.dq = propagv(mol.dq); 
     par{1}.j = propagv(mol.j);
     aux.erroz= mod( reseta(aux.erroz + propagv(zpre{1})),360);
%--------------------------------------------------------------------------
     
    for m=1:ntarg{1}
          if aux.listvirtua(targ{1}(m)) == 0; 
             aux.listvirtua(targ{1}(m)) = 1; aux.phasevirtua(targ{1}(m)) = mod(phase{1}(m),360);
          else
            aux.listvirtua(targ{1}(m)) = 0;  phase{k}(m) = aux.phasevirtua(targ{1}(m)); aux.phasevirtua(targ{k}(m)) = 0;
          end
    end
   
%--------------------------------------------------------------------------
    par{2}.j = propagv(mol.j);
    par{2}.zz = reseta(propagv(zzang + zzpos{1})); 
    par{2}.dq = propagv(mol.dq); 
    zpos{1} =  reseta(propagv(zpos{1} +zang));
 
%--------------------------------------------------------------------------
    tmax = 2/abs(mol.j(aux.atarg{1}.targ(1),aux.atarg{1}.targ(2)));
    
    x1 = fminbnd(@(x) fmin3(x,par),0,tmax,aux.opt);
    
    aa = 360*x1*mol.j(aux.atarg{1}.targ(1),aux.atarg{1}.targ(2));
   
    zzflag = 0;
    if  abs(aa) > 360;  
        
        x1 = mod(abs(aa),360)/(360 * abs(mol.j(aux.atarg{1}.targ(1),aux.atarg{1}.targ(2))));
        aux.count = aux.count +1; 
        struc{aux.count}.tip = 'rotzp';
        struc{aux.count}.desc = 0*(1:aux.N);
        struc{aux.count}.desc(aux.atarg{1}.targ(1)) = sign(aa)*180;
       
    
        aux.count = aux.count +1; 
        struc{aux.count}.tip = 'rotzp';
        struc{aux.count}.desc = 0*(1:aux.N);
        struc{aux.count}.desc(aux.atarg{1}.targ(2)) = sign(aa)*180;
        zzflag = 1;
    end
    
    x2 = 0;
    for k=1:length(aux.refoc);
         par{3}.t = x1; 
         x2 = fminbnd(@(x) fmin4(x,par),-x1/2,x1/2,aux.opt) ;
    end

   
    
%--------------------------------------------------------------------------
    if  x2 == 0;
    topt = [x1]
    aux.erroz = mod(reseta(aux.erroz + 360*par{1}.dq*topt(1)),360);
    aux.count = aux.count +1; 
	struc{aux.count}.tip = 'delay';   
	struc{aux.count}.desc = [topt(1)];
    aux.errozz = mod(reseta(par{1}.zz + par{2}.zz + 360*par{1}.j*topt(1)),720);
    else
        topt = [x1/2+x2 x1/2-x2]; 
        aux.erroz = mod(reseta(aux.erroz + 360*par{1}.dq*topt(1)),360);
        aux.count = aux.count +1; 
	    struc{aux.count}.tip = 'delay';   
	    struc{aux.count}.desc = [topt(1)];
        
        for k=1:length(aux.refoc);
    
            phasep = mod(phase{1},360);
            
            aux.count = aux.count +1; 
            struc{aux.count}.tip = 'pulse'; 
            struc{aux.count}.desc = aux.p{pul}; 
            struc{aux.count}.phase = phasep; 
            struc{aux.count}.list = aux.listspin;
            
            aux.count = aux.count +1; 
            struc{aux.count}.tip = 'rotzp';
            struc{aux.count}.desc = 0*(1:aux.N);
            struc{aux.count}.desc(targ{1}) = -aux.erroz(targ{1});
            struc{aux.count}.desc  = propagv(struc{aux.count}.desc);
            aux.erroz(targ{1}) = 0;
        end
        
        aux.erroz = mod(reseta(aux.erroz + zpos{1}  + 360*par{2}.dq*topt(2)),360);
	    aux.count = aux.count +1; 
	    struc{aux.count}.tip = 'delay';   
	    struc{aux.count}.desc = [topt(2)];
        aux.errozz = mod(reseta(par{1}.zz + par{2}.zz + 360*par{1}.j*topt(1) + 360*par{2}.j*topt(2)),720);
    end
    
%--------------------------------------------------------------------------
        if zzflag == 1
    aux.errozz(aux.atarg{1}.targ(1),aux.atarg{1}.targ(2)) = ...
        aux.errozz(aux.atarg{1}.targ(1),aux.atarg{1}.targ(2)) + sign(aa)*360;
        end
    aux.ttot = aux.ttot + sum(topt) + pw{1};    
    aux.ttotp =  aux.ttotp + pw{1};

U1 =1; U2=1;

for k=1:length(aux.atarg); 
     U1 = zzgate(aux.N,aux.atarg{k}.ang*pi/180,aux.atarg{k}.targ(1),aux.atarg{k}.targ(2))*U1;
     U2 = zzgate(aux.N,aux.errozz(aux.atarg{k}.targ(1),aux.atarg{k}.targ(2))*pi/180, ...
         aux.atarg{k}.targ(1),aux.atarg{k}.targ(2))*U2;
         aux.errozz(aux.atarg{k}.targ(1),aux.atarg{k}.targ(2)) = 0;
end
aux.Uid{aux.count2} = U1*aux.Uid{aux.count2};
aux.Usim{aux.count2} = U2*aux.Usim{aux.count2};

aux.atarg = [];
aux.atargr = [];
aux.refoc = [];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function f = fmin3(x,par)
    global mol aux
            
              nn= aux.atarg{1}.targ;
              angt = aux.atarg{1}.ang;

               ang1= reseta(par{1}.zz + 360*par{1}.j*x/2);
               ang2= reseta(par{2}.zz + 360*par{2}.j*x/2);
               ang = mod(ang1+ ang2,720);
               U = zzgate(aux.N,pi*angt/180,nn(1),nn(2)); 
               V = zzgate(aux.N,pi*ang(nn(1),nn(2))/180,nn(1),nn(2)); 
   
               f = 1 - fidelmat(U,V);
             
    end
    
function f = fmin4(x,par)
    global mol aux
               t = par{3}.t;
               ang1= reseta(par{1}.zz + 360*par{1}.j*(t/2 + x));
               ang2= reseta(par{2}.zz + 360*par{2}.j*(t/2 - x));
               ang = mod(ang1+ ang2,720);
              
               U=eye(2^aux.N); V=1;
               for k=1:length(aux.atargr{1}.targ);
                   nn= aux.atargr{1}.targ{k};
                   V = zzgate(aux.N,pi*ang(nn(1),nn(2))/180,nn(1),nn(2))*V;
               end
   
               f = 1 - fidelmat(U,V);            
 end   
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function apdelay(targ1,targ2,ang)
% Set um alvo para o delay 
global struc aux mol spectro; 
aux.flag = 1; 
aux.atarg{1}.targ = [targ1 targ2];
aux.atarg{1}.ang = ang;
end
    
%--------------------------------------------------------------------------
function aprefocus(nome,phase,targ)
%Set o alvo para arefocaliza��o
global struc aux mol  spectro; 
if aux.flag == 0; error('Refocus must be after a delay'); end
m = length(aux.refoc)+1;
aux.refoc{m}.nome = nome;
aux.refoc{m}.phase = phase;
[nn bb] = size(targ);
for k=1:nn; aux.atargr{m}.targ{k} = targ(k,:); end

end

%--------------------------------------------------------------------------

function aprz(ang,targ);
% lida com rota��es em Z
global struc aux; 

aux.count = aux.count + 1; 
struc{aux.count}.tip = 'rotzp'; 
struc{aux.count}.desc(1:aux.N) = 0; struc{aux.count}.desc(targ) = ang; 
aux.Uid{aux.count2} = multigate(aux.N,targ,rotz(ang*pi/180))*aux.Uid{aux.count2};
aux.Usim{aux.count2} = multigate(aux.N,targ,rotz(ang*pi/180))*aux.Usim{aux.count2};
struc{aux.count}.desc = propagv(struc{aux.count}.desc);

end


function  resetspin(vec); 
%atualiza listspin
global aux;  
aux.listspin(vec) = 1; 
aux.errozz=  reseta(aux.errozz); 
aux.erroz=  reseta(aux.erroz); 
end


function  unresetspin(vec); 
%atualiza listspin
global aux;  
aux.listspin(vec) = 0;
end


function  initspin(vec); 
%atualiza listspin
global aux;  
aux.listspin(vec) = 1;
aux.list0=aux.listspin;
end
%--------------------------------------------------------------------------
function  D = reseta(angin);
%zera os erros deacordo com a listspin
global aux
    n =size(angin);  D =angin;
    
    if n(1) ==1; 
        D(aux.listspin == 1) = 0; 
    else
       
       for k=1:n(2) ; for m=k+1:n(2); if aux.listspin(k) == 1 & aux.listspin(m) == 1; D(k,m) = 0 ;end; end; end
       
    end
   
end

%--------------------------------------------------------------------------
function  D = propagv(angin);
    %propaga pulsos virtuais
    global aux
    n = size(angin);  D =angin;
    
    if n(1) ==1; 
        D(aux.listvirtua == 1) = -D(aux.listvirtua == 1); 
    else

      for k=1:n(2) ;
		  
          for m=k+1:n(2) ;
			  if aux.listvirtua(k) == 1; D(k,m) = -D(k,m) ; end
              if aux.listvirtua(m) == 1; D(k,m) = -D(k,m) ;end; 
          end; 
      end 
       
    end
   
end

%--------------------------------------------------------------------------

function pul=findpul(nome,pulse);
%encontra o pulso
pul = 0;

for k=1:length(pulse); 
    if strmatch(pulse{k}.name,nome,'exact'); pul = k ;break; end
end;
    if pul == 0; error([nome ' IS NOT IN THE PULSE LIST']); end
        
end

%--------------------------------------------------------------------------

function disper(ponto);
%mostra erros
global aux

disp(['Error at ' ponto])
disp(aux.errozz)


end

function opthere;
global aux

optdelay(zeros(aux.N),zeros(1,aux.N));
   aux.flag = 0;   


end

%--------------------------------------------------------------------------

function grad(amp,temp,forma)
%lida com gradiente

global aux struc

finseq;
   
resetspin(1:aux.N);  
aux.phasevirtua(1:aux.N) = 0; 

aux.count2 = aux.count2 + 1;     
aux.Uid{aux.count2}  = eye(2^aux.N);
aux.Usim{aux.count2}  = eye(2^aux.N);

aux.count = aux.count +1;
struc{aux.count}.tip = 'grad'; 
struc{aux.count}.desc.amp = amp;
struc{aux.count}.desc.temp = temp;
if nargin == 3; struc{aux.count}.desc.forma = forma; end
end

function finseq
global aux struc

if aux.flag == 1;
   optdelay(zeros(aux.N),zeros(1,aux.N));
   aux.flag = 0;   
end

 for k=1:aux.N; 
    for m=k+1:aux.N
        aux.Usim{aux.count2} = zzgate(aux.N,aux.errozz(k,m)*pi/180,k,m)*aux.Usim{aux.count2}; 
    end
 end
 
aux.count = aux.count +1; 
struc{aux.count}.tip = 'rotzp';
struc{aux.count}.desc = propagv(-aux.erroz);
 
 
end

function strucout = writeoutput(struc);
global aux

erroz(1:aux.N) = 0;
 
m = 0;
 
for k=aux.count:-1:1
        switch struc{k}.tip
              case 'pulse'
               targ = struc{k}.desc.targ;
               struc{k}.phase  = struc{k}.phase + erroz(targ);
  			   case 'rotzp'
                  
               erroz = mod(struc{k}.desc+ erroz,360);
               case 'grad'
                 erroz(1:aux.N) = 0;
         end
  end

 erroz  = erroz.*not(aux.list0);

 for k=1:aux.count
        switch struc{k}.tip
              case 'pulse'
               targ = struc{k}.desc.targ;
               struc{k}.phase  = struc{k}.phase - erroz(targ);
               m = m+1; strucout{m} = struc{k};
               case 'delay'
                 m = m+1; strucout{m} = struc{k};
               case 'grad'
                 m = m+1; strucout{m} = struc{k};
                 erroz(1:aux.N) = 0;
         end
  end

      
strucout{m+1}.rz = mod(erroz,360);
strucout{m+1}.virtua = aux.listvirtua;
strucout{m+1}.phasevirtua = aux.phasevirtua;
strucout{m+1}.Uid = aux.Uid;

end

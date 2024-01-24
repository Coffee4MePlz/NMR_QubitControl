
function  strucout = compilaSeqB(seq,pulse,seq2);
%strucout = compilaSeqB(seq,pulse);
%Compile a sequence "seq" 
%
clear  aux struc;
global mol spectro aux struc;

aux.opt = optimset('TolFun',1e-9,'MaxFunEvals',5e10,'TolX',1e-9);
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
%fun��o lida com um pulso, 

global aux struc;
      
   %Procura o pulso  
    pul=findpul(nome,aux.p);   
    
    %Ler os par�metros do pulse 
    targ = aux.p{pul}.targ;% Spin Alvo
    ang = aux.p{pul}.ang;%Angulo aplicado 
    pw = aux.p{pul}.pw;% Tempo do pulso
    zzpos = aux.p{pul}.zzp; %Pos-ZZ erro
    zzpre = aux.p{pul}.zza;%Pre-ZZ erro 
    zpre = aux.p{pul}.za;%pre-Z erro
    zpos = aux.p{pul}.zp; %pos-Z erro
    ntarg = length(targ);
    
%--------------------------------------------------------------------------
%Se n�o hover um delay antes do pulso: apenas atualiza os erros com o pre-error 
%do pulse. Sen�o otimiza o delay de acordo com o foi pedido 
if aux.flag == 0;
    aux.errozz = reseta(aux.errozz + propagv(zzpre));
    aux.erroz= mod( reseta(aux.erroz + propagv(zpre)),360);
	%atualiza o propagador da simula��o at� antes do pulso
else
   optdelay(zzpre,zpre);
   aux.flag = 0;
      
end


aux.count = aux.count +1; 
struc{aux.count}.tip = 'rotzp';
struc{aux.count}.desc = 0*(1:aux.N);
struc{aux.count}.desc(targ) = -aux.erroz(targ);
struc{aux.count}.desc  = propagv(struc{aux.count}.desc);
aux.erroz(targ) = 0;



%Zera erros referente ao targ antes do pulso
for m=1:ntarg
   
for k=1:aux.N; aux.Usim{aux.count2} = zzgate(aux.N,(pi/180)*(aux.errozz(targ(m),k) +...
        aux.errozz(k,targ(m))),targ(m),k)*aux.Usim{aux.count2}; end
    
aux.errozz(targ(m),:) = 0; aux.errozz(:,targ(m)) = 0;

aux.listspin(targ(m)) = 0;
end


if nargin == 3; aux.listspin(optin) = 1; end

%Se houver pulso virtual no targ, combina com p pulso

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

%escreve um novo evento
for m=1:ntarg
aux.Usim{aux.count2} = multigate(aux.N,targ(m),rotphi(aux.p{pul}.ang(m)*pi/180,phase(m)*pi/180))*aux.Usim{aux.count2};
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


%atualiza erros, propaga Z para o fim e propaga pulsos virtuais para o fim 
aux.errozz = reseta(aux.errozz + propagv(zzpos));

aux.erroz= mod( reseta(aux.erroz + propagv(zpos)),360);

%atualiza o tempo total
aux.ttot = aux.ttot + pw;  
aux.ttotp = aux.ttotp  + pw;
%atualza os propagadores
for m=1:ntarg
aux.Uid{aux.count2} = multigate(aux.N,targ(m),rotphi(aux.p{pul}.ang(m)*pi/180,phase(m)*pi/180))*aux.Uid{aux.count2};
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function optdelay(zzang,zang)
%funcao otimiza delays deacordo com o que foi pedido 

global mol spectro aux struc;

if isempty(aux.refoc)
    %propaga errors ZZ e pulsos virtuais at� o fim
	
	par.zz = reseta(aux.errozz + propagv(zzang ));
    
    par.j = propagv(mol.j);
	
	for m=1:length(aux.atarg)
    
		nn = aux.atarg{m}.targ;
        aa = mod(aux.atarg{m}.ang - par.zz(nn(1),nn(2)),720);
	
	    tt1(m) = aa/(360*par.j(nn(1),nn(2)));

	    if tt1(m) < 0; aa = mod(aa,-720); tt1(m) = aa/(360*par.j(nn(1),nn(2))); end  
	
	end
	topt = max(tt1);
	
	%atualiza erros, propaga Z para o fim e propaga pulsos virtuais para o fim 
    aux.erroz = mod(reseta(aux.erroz + propagv(zang ) + 360*propagv(mol.dq)*topt ),360);
    
    aux.errozz = reseta(par.zz + 360*par.j*topt);
    
    %escreve um novo evento
    aux.count = aux.count +1; struc{aux.count}.tip = 'delay';   struc{aux.count}.desc = [topt];
    
    %atualiza o tempo total
    aux.ttot = aux.ttot + sum(topt);  
    
else
	%par�metros do pulso de refocaloza��o
    pul=findpul(aux.refoc.nome,aux.p);   
    targ = aux.p{pul}.targ; 
    ang = aux.p{pul}.ang;
    pw = aux.p{pul}.pw;     
    zzpos = aux.p{pul}.zzp;     
    zzpre = aux.p{pul}.zza;
    zpre = aux.p{pul}.za;   
    zpos = aux.p{pul}.zp;  
    phase = aux.refoc.phase;
    ntarg = length(targ);
    
	%acho tempo total do delay
	
	%clear par;
    %propaga Z e pulsos virtuais at� o pulso de refocaliza��o 
    par{1}.zz = reseta(aux.errozz + propagv(zzpre)); 
    par{1}.j = propagv(mol.j); 
    par{1}.dq = propagv(mol.dq); 
    aux.erroz= mod( reseta(aux.erroz + propagv(zpre)),360);
    %Se houver pulso virtual do targ combina o pulso de refocaliza��o ou
    %cria um novo
    
    for m=1:ntarg
    
    if aux.listvirtua(targ(m)) == 0; 
        aux.listvirtua(targ(m)) = 1; aux.phasevirtua(targ(m)) = mod(phase(m),360);
    else
      
        aux.listvirtua(targ(m)) = 0;  phase(m) = aux.phasevirtua(targ(m)); aux.phasevirtua(targ(m)) = 0;
    end
    
    end
    
    %efeito do pulso de refocaliza��o no ZZ e no Z
    par{2}.j = propagv(mol.j);
    par{2}.zz = reseta(propagv(zzang + zzpos)); 
    par{2}.dq = propagv(mol.dq); 
    zpos =  reseta(propagv(zpos +zang));
 
    %encontra a posicao do refoc
	for m=1:length(aux.atarg)
    
		nn = aux.atarg{m}.targ;
        aa = mod(aux.atarg{m}.ang - par{1}.zz(nn(1),nn(2)) ...
		- par{2}.zz(nn(1),nn(2)),720);
	 
	    tt1(m) = aa/(360*par{1}.j(nn(1),nn(2)));

	    if tt1(m) < 0; aa = mod(aa,-720); tt1(m) = aa/(360*par{1}.j(nn(1),nn(2))); end  
	
    end
    
	topt1 = max(tt1);par{end+1} = topt1; 
	
	topt2 = fminbnd(@(x) fmin2(x,par),-topt1/2,topt1/2,aux.opt);  
   
	
    topt = [topt1/2-topt2 topt1/2+topt2];
     %erros devido ao promeiro delay
    aux.erroz = mod(reseta(aux.erroz + 360*par{1}.dq*topt(1)),360);
    
    %escreve o delay
    aux.count = aux.count +1; 
	struc{aux.count}.tip = 'delay';   
	struc{aux.count}.desc = [topt(1)];
    
	
    %escrve o pulso
    phasep = mod(phase,360);
	
   
    aux.count = aux.count +1; 
    struc{aux.count}.tip = 'pulse'; 
    struc{aux.count}.desc = aux.p{pul}; 
    struc{aux.count}.phase = phasep; 
    struc{aux.count}.list = aux.listspin;
 
    aux.count = aux.count +1; 
    struc{aux.count}.tip = 'rotzp';
    struc{aux.count}.desc = 0*(1:aux.N);
    struc{aux.count}.desc(targ) = -aux.erroz(targ);
    struc{aux.count}.desc  = propagv(struc{aux.count}.desc);
    aux.erroz(targ) = 0;
    
   
	%errors devido ao segundo delay
    aux.erroz = mod(reseta(aux.erroz + zpos  + 360*par{2}.dq*topt(2)),360);
	
    %escreve o segundo delay
    aux.count = aux.count +1; 
	struc{aux.count}.tip = 'delay';   
	struc{aux.count}.desc = [topt(2)];
 
	
    %atualiza error ZZ
    aux.errozz = reseta(par{1}.zz + par{2}.zz +  360*par{1}.j*topt(1) + 360*par{2}.j*topt(2));
    
    %atualiza o tempo total
    aux.ttot = aux.ttot + sum(topt) + pw;    
    aux.ttotp =  aux.ttotp + pw;
end

%atualiza os propagadores
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

    function f = fmin2(x,par)
    %encontra um delay um delay         
    global mol aux

  t1 = par{end};
	
	
    ang= reseta(par{1}.zz + 360*par{1}.j*(t1/2-x));
	ang2= reseta(par{2}.zz + 360*par{2}.j*(t1/2+x)); 
	
    ang = mod(ang + ang2,720);
	
	for m=1:length(aux.atarg)
    mm= aux.atarg{m}.targ;
    ang(mm(1),mm(2)) = aux.atarg{m}.ang;
	end
	

    U =1;
    for k=1:aux.N; 
		for n=k+1:aux.N; 
			U = zzgate(aux.N,pi*ang(k,n)/180,k,n)*U; 
		end;
    end
    
   
    bb = length(aux.atargr{1}.targ);
    
    V = U;
    
	nn = aux.atargr{1}.targ;
	
	for k=1:bb
    
		
	V = zzgate(aux.N,-pi*ang(nn{k}(1),nn{k}(2))/180,nn{k}(1),nn{k}(2))*V;
	
    end
    
    f = 1 - fidelmat(U,V);
             
 end   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function apdelay(targ1,targ2,ang)
% Set um alvo para o delay 
global struc aux mol spectro; 

aux.flag = 1; 
m=length(aux.atarg);
aux.atarg{m+1}.targ = [targ1 targ2];
aux.atarg{m+1}.ang = ang;

end
    
%--------------------------------------------------------------------------
function aprefocus(nome,phase,targ)
%Set o alvo para arefocaliza��o
global struc aux mol  spectro; 

if aux.flag == 0; error('Refocus must be after a delay'); end

aux.refoc.nome = nome;
aux.refoc.phase = phase;
[nn bb] = size(targ);
for k=1:nn
aux.atargr{1}.targ{k} = targ(k,:);
end

end

%--------------------------------------------------------------------------
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

    
    %----------------------------------------------------------------------
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

function grad(amp,temp)
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


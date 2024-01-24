function roh = tomography(pul,struc,nq,tomospin,spinobs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Tomography for N spins (isisng coupling)                                                                   %          
%Alexandre M . de Souza                                                                 %
%function roh = tomography(pul,struc,nq,tomospin);                                      %
%pul = preparatory pulses EX: {+II,+IX,+IY,+XX}                                    %
%struc = structure with fitted data struc{k,m} (spin k, experiment m)         %
%nq = nunber of qubtis                                                                    %
%tomospin = qubits not tomographed (optional)                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BASE PARA EXPANDIR A MATRIZ DENSIDADE                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 4
      
       tomospin = sort(tomospin);
       aux =1:nq; 
       aux(tomospin) = [];
       [subbasiss] = subsigbasis(nq,aux);
else 
       [subbasiss] = subsigbasis(nq) ;
       tomospin=1:nq;
end

[basisn] = sigbasis(length(tomospin));
[subbasiss] = subsigbasis(nq) ;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transitionstomo ENCONTRA QUE ELEMENTOS DA BASE OBSERVADOS NO ESPECTRO DE CADA SPIN     % 
%opR - > Operadores da Parte Real  e opI - > Operadores da parte Imaginária              % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[opR{1} opI{1}] = transitionstomo(nq,spinobs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENCONTRA QUE ELEMENTOS DA BASE OBSERVADOS NO ESPECTRO DE CADA SPIN PARA CADA PULSO     % 
% listobR e listobI - > Listas de operadores,                                            %
%ssR e ssI -> lista contendo que linhas devem ser somadas em caso de tomografia parcial  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ll =1:length(pul)
  
        % Find the opearators that can be observed after the reading pulses
          opRpul{1,ll}  = transitionspul3(opR{1},pul{ll},nq,1);
          opIpul{1,ll}  = transitionspul3(opI{1},pul{ll},nq,0);
         
          [listobR{1,ll} ssR{1,ll}] = subcompare2sigb(opRpul{1,ll},subbasiss,nq);
          [listobI{1,ll} ssI{1,ll}] = subcompare2sigb(opIpul{1,ll},subbasiss,nq);
       

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RECONTRUCAO DA MATRIZ DENSIDADE                                                              % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ll=1:length(listobR);

   % for k=1:length(tomospin)
     k = 1;   
       %Pega as amplitudes ajustadas      
       [specr aaa] =  subauxxx2(struc{k,ll},ssR{k,ll});
       [aaa speci] =  subauxxx2(struc{k,ll},ssI{k,ll});
        
       %Encontra os coeficientes de cada elemnto da base apartir das amplitudes
       matelemenR{k,ll} = solveop(specr,listobR{k,ll},ssR{k,ll},nq);
       matelemenI{k,ll} = solveop(speci,listobI{k,ll},ssI{k,ll},nq);

    %end
    
end

%Combina todos os coeficientes para formar a matriz
roh=combineall(matelemenR,matelemenI,listobR,listobI,basisn);

roh=ztrace(multigate(nq,spinobs,xgate)*roh,spinobs);

disp(roh)

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [eer eei] = subauxxx2(ss,sss);
%funcao para pegar as amplitudes ajustadas

for k=1:length(sss)

    for m=1:length(sss{k})
       auxr(m) =  sum(ss(2*sss{k}(m)-1));
       auxi(m) =  sum(ss(2*sss{k}(m)));
    end
       
   % eei(k) = -sum(auxi);     eer(k) = sum(auxr); %Bruker
   %  eei(k) = sum(auxi);        eer(k) = sum(auxr); X-iY X+Y
  eei(k) = sum(auxi);        eer(k) =  - sum(auxr) ;%Varian  X+iY X+Y
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ee = solveop(spec,listob,ss,nq);
% function ee = solveop(spec,listob);
% Resolve a equações para encontrar os coeficientes assocoados aos
% operadores 

for k=1:2^(nq-1)

    for m=1:2^(nq-1)
                AA(k,m) = listob{ss{k}(1)}(2*m);  
    end
end
ee = (AA\spec');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [roh]=combineall(matelemenR,matelemenI,listobR,listobI,basisn);
% function [roh]=combineall(matelemenR,matelemenI,listobR,listobI,basisn);
% combina todos os operadores observados e monta a matriz densidade

[aa bb]=size(listobR);

vec(1:length(basisn)) = 0;

nq = log2(length(basisn))/2;
for k=1:length(basisn); aux{k} = 0 ; end
for k=1:aa
   
    for m=1:bb
         
         for kk=1:1:length(listobR{k,m}{1})/2
             vec(listobR{k,m}{1}(2*kk-1)) = vec(listobR{k,m}{1}(2*kk-1)) +1;
             ind = vec(listobR{k,m}{1}(2*kk-1)); 
             %aux{listobR{k,m}{1}(2*kk-1)}(ind) = matelemenR{k,m}(kk); 
              aux{listobR{k,m}{1}(2*kk-1)}(ind) =- matelemenR{k,m}(kk); %varian
          end
         
         for kk=1:1:length(listobI{k,m}{1})/2
             vec(listobI{k,m}{1}(2*kk-1)) = vec(listobI{k,m}{1}(2*kk-1)) +1;
             ind = vec(listobI{k,m}{1}(2*kk-1)); 
            % aux{listobI{k,m}{1}(2*kk-1)}(ind) = matelemenI{k,m}(kk); 
              aux{listobI{k,m}{1}(2*kk-1)}(ind) =- matelemenI{k,m}(kk); %varian
         end
       
         aux{1} = 0/(2^nq);
        
    end
end

roh =  aux{1} *basisn{1};

for k=2:length(basisn);
    roh = roh + mean(aux{k}) *basisn{k};
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bb] = sigbasis(nq);
% [bb] = sigbasis(nq);
% bb -> operadores de base para expandir rho
% nq - > número de qbits

op = {eye(2), [0 1;1 0], [0 -i;i 0],[1 0; 0 -1]} ;

for k=0:4^nq-1;
    
    aux =dec2base(k,4,nq);
   
    bb{k+1} =1;
       
    for m=1:nq; 
  
        aux2{m} = multigate(nq, m, op{str2num(aux(m)) + 1});     
        
        bb{k+1} = (bb{k+1}*aux2{m});
           
    end
         
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bb2] = subsigbasis(nq,sub);
% [bb bb2] = sigbasis(nq);
% bb2 -> operadores simbolicos
% nq - > número de qbits
% sub- > subsistema para ser descartado 

op2 = {'I','X','Y','Z'};

if nargin == 1;  nqs = nq; sub = 0; else; nqs = nq-length(sub); end

for k=0:4^(nqs)-1;
    
    aux =dec2base(k,4,nqs);
   
    bb2{k+1} =[];
       
    mm = 1;
    
    for m=1:nq; 
  
        if any(m*ones(1,length(sub)) == sub);
            OPs = op2{1};
        else
            OPs = op2{str2num(aux(mm)) + 1};
            mm = mm+1;
        end
              
        bb2{k+1} = [bb2{k+1} OPs];
    
    end
       
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [opR opI] = transitionstomo(N,sp)
% [tran ff T2 opR opI] = transitions(N,sp)
% Calculate the operators associated to each line of the spectrum 
% opR{k}  - > Operators associated to the line k in Real part 
% opI  - > Operators associated to the line k in Real part 
% N - > Number of spins
% sp -> Spin to caculate the the operators   
% 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transitions
% st1 is all possible states
% st2 is equal to st1 with the spin (sp) inverted
mm =1;

for k=0:2^N-1; 
    
    st = dec2bin(k,N);  
    
    if st(sp) == '0'; 
        st1{mm} = st;   
        st2{mm} = st; 
        st2{mm}(sp) = '1'; 
        mm = mm +1;
    end
   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%frequency
%If the neighour spin is in state 0 the frequecy in shifted by -J/2
%If the neighour spin is in state 1 the frequecy in shifted by +J/2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Operators Using rules of PHYSICAL REVIEW A 69, 052302 (2004)
%

 
 for k =1:length(st1)
     
     for m=1:N;
         
        if m == sp 
            
            aux2{k}(m) = 'X';  auy2{k}(m) = 'Y';
            
        else
            
            if st1{k}(m) == '0'; 
                aux2{k}(m) = 'I';  auy2{k}(m) = 'I';
            else 
                aux2{k}(m) = 'Z'; auy2{k}(m) = 'Z';
            end
            
        end
            
     end
     
    %---------------------------------------------------------------------

 end

 HH = hadamard(length(st1));
 
 for k =1:length(st1)
     
    for m = 1:length(st1)
     
        if HH(k,m) == 1; 
            
            sinalx{k}(m) = '+'; sinaly{k}(m) = '-'; 
           
        else 
            sinalx{k}(m) = '-'; sinaly{k}(m) = '+'; 
           
        end
        
 end
     
 end

for k =1:length(st1)
    
    opR{k} = [sinalx{k}(1) aux2{1}];
    
    opI{k} = [sinaly{k}(1) auy2{1}];

    for m =2:length(st1)
        
        opR{k} = [opR{k} sinalx{k}(m)  aux2{m}];
        
        opI{k} = [opI{k} sinaly{k}(m)  auy2{m}];
    end

end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
function tran2 = transitionspul3(tran,pul,nq,opt)
% tran2 = transitionspul3(tran,pul,nq)
% trans2 -> Operators associate to each line after a reading pulse 
% trans - >  Operators associate to each line before a reading pulse 
% pul - > pulse
% nq - > number os qubits 

HH = hadamard(length(tran));

for k=1:length(tran)
    
    for n=0:(length(tran) -1)
     
    signn = 1; 
    if pul(1) == '-'; signn = not(signn); end 

    for m=2:nq+1

        if tran{k}(m+n*(nq+1)) ==  'I'; tran2{k}(m+n*(nq+1)) = 'I'; end
        
        
         if tran{k}(m+n*(nq+1)) ==  'X'; 
            if pul(m) == 'X'; tran2{k}(m+n*(nq+1)) = 'X'; end
            if pul(m) == 'Y'; tran2{k}(m+n*(nq+1)) = 'Z'; end
            if pul(m) == 'I'; tran2{k}(m+n*(nq+1)) = 'X'; end
        end
        
         if tran{k}(m+n*(nq+1)) ==  'Y'; 
            if pul(m) == 'Y'; tran2{k}(m+n*(nq+1)) = 'Y'; end
            if pul(m) == 'X'; tran2{k}(m+n*(nq+1)) = 'Z';  signn = not(signn); end
            if pul(m) == 'I'; tran2{k}(m+n*(nq+1)) = 'Y'; end
         end
        
          if tran{k}(m+n*(nq+1)) ==  'Z'; 
            if pul(m) == 'X'; tran2{k}(m+n*(nq+1)) = 'Y'; end
            if pul(m) == 'Y'; tran2{k}(m+n*(nq+1)) = 'X'; signn = not(signn); end
            if pul(m) == 'I'; tran2{k}(m+n*(nq+1)) = 'Z'; end
        end
        
    end
    
    if opt == 1;
    if HH(k,n+1) == -1; signn = not(signn); end
    else 
    if HH(k,n+1) == +1; signn = not(signn); end  
    end
   
    if signn == 1; tran2{k}(1+n*(nq+1)) = '+'; else tran2{k}(1+n*(nq+1)) = '-'; end
    end   
    
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [listob ss] = subcompare2sigb(op,sigb,nq);
% [listob2 ss] = subcompare2sigb(op,sigb);
%  
% ss - >  lista contendo as linhas que devem ser adicionadas
% listob2{k} - >  list with observables  
% sigb - >  Base of subsystem to the tomographed
% op  - >   Operators associated to each line
% nq - > nnumber of spins

nconfig = [];

for k=1:length(op)
    
   n0 = 2;
   nf = n0+ nq-1;
   for m = 1:length(op)
   aux{k,m} = op{k}(n0:nf);
   sigaux{k,m} =  op{k}(n0-1);
   n0 = nf + 2;
   nf = n0 + nq-1;
   end
   
   kk = 1;
   
   for m=1:length(op)
            for n=1:length(sigb)
                    if aux{k,m} == sigb{n};
                      listob{k}(kk) = n;
                      kk=kk+1;
                      if sigaux{k,m} == '+'; listob{k}(kk) = +1; else;  listob{k}(kk) = -1; end
                      kk=kk+1;
                      break;
                    end
            end
   end
   
   flag = 0;
   [aa bb] = size(nconfig);
   
   for g=1:aa
   if listob{k}(2:2:end) ==  nconfig(g,:); ss{g} = [ss{g}(1,:) k ]; flag = 1; break; end
   end
   
   if flag == 0; 
       nconfig(aa+1,:) = listob{k}(2:2:end); 
       ss{aa+1} = [k]; 
   end

end

end
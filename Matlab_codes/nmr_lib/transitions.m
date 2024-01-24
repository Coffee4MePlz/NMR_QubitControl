function [freq tran ff opR opI] = transitions(mol,sp)
% [tran ff opR opI] = transitions(Def,sp)
% Find the frequncies and transtions of the spin "sp" of the molecule "mol"
%  tran -> transitions
%  freq -> frequency
%  opR, opI - > product operators assossiated to each transition 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transitions
% st1 is all possible states
% st2 is equal to st1 with the spin (sp) inverted

N = mol.nspin;  

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

for k=1:length(st1); tran{k} = [st1{k} ' <--> ' st2{k}]; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%frequency
%If the neighour is in the state 0 the frequecy in shifted by -J/2
%If the neighour is in the state 1 the frequecy in shifted by +J/2

for k=1:length(st1);
    
    freq{k} = [ 'W(' mat2str(sp) ')'];
  
    for m=1:N;
        if m ~= sp 
            
            %if st1{k}(m) == '0'; sinal(m) = '-';else; sinal(m) = '+'; end %Broker
            if st1{k}(m) == '0'; sinal(m) = '+';else; sinal(m) = '-'; end  %VARIAN
            jj = mat2str(sort([sp m]));
            
            freq{k} = [ freq{k}  sinal(m) 'J' jj '/2'  ];
            
        end
        
    end
    
end

%Read molucule file to find numerical values to freqs
   
T2 = mol.T2(sp);

JJ = mol.j;

wmega = mol.dq(sp);

for m=1:N; 
    k=(m+1):N; 
    JJ(k,m) = mol.j(m,k); 
    if JJ(k,m) > 1; sigj(k,m) = 1; else sigj(k,m) = -1; end;
end;

aux = 1:N; aux(aux==sp) = [];

for k=1:length(freq);
   
    ff(k) = wmega;
    
    sigmm = freq{k}(5:9:end);

    for m=1:N-1
 
        
              if sigmm(m) == '+'; 
                          ff(k) = ff(k) + JJ(sp,aux(m))/2; 
              else;  
                          ff(k) = ff(k) - JJ(sp,aux(m))/2;
              end     

    end
    
end

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
     

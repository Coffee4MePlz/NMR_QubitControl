function struc = bruker2mat2(cca,v,nump,freg);
%to read bruker files  
% C:\Users\John\Documents\MATLAB\dadosfid\Chloroform'
% cca = nome da pasta
% v = numero da pasta do fid
% nump = mumero do experimento (pasta dentro da pasta pdata

%  cca = cca


for k=1:length(v) 

    cc = [cca '\' num2str(v(k)) '\'];
            
%Reading ACQUS------------------------------------------------------------- 
          
    aux1 = fopen([cc 'acqus']);
            
 if aux1 == -1; 
    
     error('files "acqus" not found in current path  I can�t proceed')
     
 else 
               
     straux2 = [];
     
               
     for  m=1:405
                    
         straux = fgetl(aux1); 
                    
         straux([findstr(straux, '#') findstr(straux, '$')]) = '';
                    
         straux2 = str2mat(straux2,deblank(straux));
         
               
     end
     
  %Solvent
  sstring  =  straux2(strmatch('SOLVENT'  , straux2),:);
  struc.sim.par{k}.SOLVENT = sstring(10:20);   
  
   %Nuc
   sstring  =  straux2(strmatch('NUC1'  , straux2),:);
   struc.sim.par{k}.nuc = sstring(7:11);   
     
 %Modo de Aquisicao    
 eval([straux2(strmatch('AQ_mod='  , straux2),:) ';']);
            
 struc.sim.par{k}.aq_mod = AQ_mod;     
      

 %ASCII FORMAT for FID
            
 eval([straux2(strmatch('BYTORDA='  , straux2),:) ';']); 

 struc.sim.par{k}.BYTORDA = BYTORDA;
 
 
 %FID scaling factor 
              
 eval([straux2(strmatch('NC='  , straux2),:) ';']); 
           
 struc.sim.par{k}.NC = NC;

 
 %Number of points in FID
            
 eval([straux2(strmatch('TD='  , straux2),:) ';']);
            
 struc.sim.par{k}.TD = TD/2; 
       
            
 %Spectral Window (Hz)
            
 eval([straux2(strmatch('SW_h='  , straux2),:) ';']);
            
 struc.sim.par{k}.sw = SW_h;
       
%             
 %Temperature (Kelvin)
            
 eval([straux2(strmatch('TE='  , straux2),:) ';']);
            
 struc.sim.par{k}.Temp = TE; 
       
            
 %Delay beetween last pulse and begining of aquisition (microseconds)
            
 eval([straux2(strmatch('DE='  , straux2),:) ';']);
            
 struc.sim.par{k}.DE = DE; 
       
 %DECIM
            
 eval([straux2(strmatch('DECIM='  , straux2),:) ';']);
            
 struc.sim.par{k}.DECIM = DECIM; 
  
 %DSPFVS
            
 eval([straux2(strmatch('DSPFVS='  , straux2),:) ';']);
            
 struc.sim.par{k}.DSPFVS = DSPFVS; 
            
 end
            
        fclose(aux1);
       
            
%Reading Fid---------------------------------------------------------------
            
aux1 = fopen([cc 'fid']);
            
if aux1 == -1; 
               
    error('files "fid" not found in current path  I can�t proceed')
            
else 
                
    if BYTORDA == 0; ascciitype = 'l';end; 
    if BYTORDA == 1; ascciitype = 'b';end;
                    
    fid  =  fread(aux1,'int32',ascciitype)';
                    
    fidr = detrend(fid(1:2:length(fid)));          
    
    fidi = detrend(fid(2:2:length(fid)));
                    
    struc.sim.fid{k} = (fidr + i*fidi)*2^(NC); 
    
end
    
            fclose(aux1);
            
           
%Reading Spectrum----------------------------------------------------------
            
aux1 = fopen([cc 'pdata\' num2str(nump) '\1r']); 
            
aux2 = fopen([cc 'pdata\' num2str(nump) '\1i']);
            
if aux1 ~= -1 & aux2 ~= -1; 
               
%Reading Procs ------------------------------------------------------------
                
aux3 = fopen([cc 'pdata\' num2str(nump) '\procs']);
                
if aux3 == -1; 
                   
    error('files "procs" not found in current path  I can�t proceed')
                
else 
                  
    straux2 = [];
                  
    for  m=1:109
                    
        straux = fgetl(aux3); 
                    
        straux([findstr(straux, '#') findstr(straux, '$')]) = '';
                    
        straux2 = str2mat(straux2,deblank(straux));
                  
    end
    
    %FT_mod
    
    eval([straux2(strmatch('FT_mod='  , straux2),:) ';']);
               
    struc.sim.par{k}.FT_mod = FT_mod;   
    
    %BC_mod
    
    eval([straux2(strmatch('BC_mod='  , straux2),:) ';']);
               
    struc.sim.par{k}.BC_mod = BC_mod;
    
    %PH_mod
    
    eval([straux2(strmatch('PH_mod='  , straux2),:) ';']);
               
    struc.sim.par{k}.PH_mod = PH_mod;  
    
    %WDW
    
    eval([straux2(strmatch('WDW='  , straux2),:) ';']);
               
    struc.sim.par{k}.WDW = WDW;  
%     
    %ASCCII FORMAT for Spectrum
               
    eval([straux2(strmatch('BYTORDP='  , straux2),:) ';']); 
    
    struc.sim.par{k}.BYTORDP = BYTORDP;
    
    %Spectrum scaling factor 
               
    eval([straux2(strmatch('NC_proc='  , straux2),:) ';']); 
              
    struc.sim.par{k}.NC_proc = NC_proc;
    
    %Number of Points in Spectrum
               
    eval([straux2(strmatch('SI='  , straux2),:) ';']);
               
    struc.sim.par{k}.SI = SI;   
    
    %Number of Points in Spectrum
               
    eval([straux2(strmatch('LB='  , straux2),:) ';']);
               
    struc.sim.par{k}.LB = LB;   
    
    %Signal to Noise Ratio
               
    eval([straux2(strmatch('SINO='  , straux2),:) ';']); 
               
    struc.sim.par{k}.signois = SINO;  
              
end
                
                
    if BYTORDP == 0; ascciitype = 'l'; end; if BYTORDP == 1; ascciitype = 'b';end;
                
    er  =  detrend(fread(aux1,'int32',ascciitype))';    
    
    ei  =  detrend(fread(aux2,'int32',ascciitype))';
           
    
    if nargin == 3;
        dt = 1/SW_h;
        nesp = length(er);
        x1 = round(nesp*dt*(freg(1) + SW_h/2));
        x2 = round(nesp*dt*(freg(2) + SW_h/2));
        fr = (1/dt)*(0:nesp-1)/nesp - SW_h/2;
        struc.sim.esp{k} = ((er(x1:x2) + i*ei(x1:x2))*2^(NC_proc)); 
        SWW = (freg(2) - freg(1));
    else
        struc.sim.esp{k} = ((er + i*ei)*2^(NC_proc)); 
        SWW = SW_h;
    end  
            
    fclose(aux1); 
            
    fclose(aux2);
            
    fclose(aux3);
            
end
    struc.sw = SWW;   
    struc.nuc = struc.sim.par{1}.nuc;
    
end


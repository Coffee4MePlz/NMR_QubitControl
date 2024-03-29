
function seq2spec(nome,seq,pul);

    global spectro mol 
    
	 fclose all;
	 ddir = '/home/quantum/vnmrsys/psglib/'; 
	 ddir2 = '/home/quantum/vnmrsys/shapelibcraph/'; 
     ddir3 = '/home/quantum/vnmrsys/maclib/'; 
	
     saida = fopen( [ddir nome '.c'],'w');
     saida2 = fopen( [ddir3 nome 'mac'],'w');
     
     fprintf(saida,'#include <standard.h> \n');
     fprintf(saida,' \n \n \n'); 
     fprintf(saida,'pulsesequence() \n');
     fprintf(saida,'{ \n'); 
     fprintf(saida,' \n');
     fprintf(saida,'double tpwr, dpwr, tof, dof, rof1, rof2, pw, p1, d1, gt1, gzlvl1; \n'); 
     fprintf(saida,'int  table1[1] = {0};\n'); 
     for k=1:length(pul); 
          for m = 1:length(pul{k}.chan)
          fprintf(saida,['double pwrf' pul{k}.name 'Ch' mat2str(pul{k}.chan{m})  ';\n']); 
          end
     end
     
	 fprintf(saida,'\n');
     fprintf(saida,'pw = getval("pw"); \n'); 
     fprintf(saida,'p1 = getval("p1"); \n'); 
     fprintf(saida,'d1 = getval("d1"); \n'); 
	 fprintf(saida,'rof1 = getval("rof1"); \n'); 
	 fprintf(saida,'rof2 = getval("rof2"); \n'); 
	 fprintf(saida,'tpwr = getval("tpwr"); \n'); 
	 fprintf(saida,'dpwr = getval("dpwr"); \n'); 
	 fprintf(saida,'tof = getval("tof"); \n'); 
	 fprintf(saida,'dof = getval("dof"); \n'); 
	 
	 for k=1:length(pul); 
         for m = 1:length(pul{k}.chan)
         sstring = ['pwrf' pul{k}.name 'Ch' mat2str(pul{k}.chan{m})];
		 fprintf(saida,[sstring '= getval("' sstring '");\n']);
         fprintf(saida2,['create(''' sstring ''',''integer'')\n']);
         fprintf(saida2,['setlimit(''' sstring ''',4095,0,1)\n']);
         end
	 end

      fprintf(saida,'\n');

	  phr = round(mod(seq{end}.rz(spectro.obs),360)/0.25);

	  fprintf(saida,'rcvrstepsize(0.25); \n');
	  fprintf(saida,'initval(%4.1f,v1); \n',phr);
	  fprintf(saida,'rcvrphase(v1); \n');
      fprintf(saida,'settable(t1,1,table1); \n');
      fprintf(saida,'setreceiver(t1); \n');
      fprintf(saida,'obspower(tpwr); \n');
      fprintf(saida,'decpower(dpwr); \n');
      fprintf(saida,'obsoffset(tof);  \n');
      fprintf(saida,'decoffset(dof);  \n \n');  
      fprintf(saida,'gt1  = 0.002;\n');   
      fprintf(saida,'gzlvl1 = 1234.0;\n');   

      fprintf(saida,'\n');
       fprintf(saida,'lk_sample(); \n');
       fprintf(saida,'delay(d1-5); \n');
       fprintf(saida,'lk_hold(); \n');
 	  fprintf(saida,'delay(5); \n');
%      fprintf(saida,'delay(d1); \n');
	  fprintf(saida,'\n');

   
for k=1:length(seq)-1
    
    switch seq{k}.tip
               case 'pulse'
                         
                         nchan  = length(seq{k}.desc.chan);
                         
                         for m=1:nchan
                         
                         n(m)=seq{k}.desc.chan{m};
                         
                         nn{m} = seq{k}.desc.formaspec{m};
                         
                          
						 aa  =  strfind(nn{m},'.RF');
						 
						 nn{m}(aa:end) = [];
						 end
                         
                          if nchan == 1;
                            if n == 1;
                             fprintf(saida,['obspwrf(pwrf' seq{k}.desc.name 'Ch1);\n']); 
                             fprintf(saida,['shapedpulse("' nn{1} '",%7.6f,zero,rof1,rof2);\n'],seq{k}.desc.pw);
                            
                            else
                            fprintf(saida,['decpwrf(pwrf' seq{k}.desc.name 'Ch2);\n']);   
                            fprintf(saida,['decshapedpulse("' nn{1} '",%7.6f,zero,rof1,rof2);\n'],seq{k}.desc.pw);    
                                
                            end
                           
                          else
                          if n(1) == 1; n1 =1; n2 = 2; else n1 = 2; n2 =1; end;
                          fprintf(saida,['obspwrf(pwrf' seq{k}.desc.name 'Ch1);\n']);   
                          fprintf(saida,['decpwrf(pwrf' seq{k}.desc.name 'Ch2);\n']);   
                          fprintf(saida,['simshapedpulse("' nn{n1} ',' nn{n2} '",%7.6f,zero,zero,rof1,rof2);\n'],seq{k}.desc.pw,seq{k}.desc.pw);
                          end
                          
               case 'delay'
              
               fprintf(saida,['delay(%8.7f);\n'],seq{k}.desc);
               
               case 'grad'
				  if isfield(seq{k}.desc,'forma')
				  fprintf(saida,['shapedgradient("' seq{k}.desc.forma '",' ...
					  mat2str(seq{k}.desc.temp) ',' ...
					  mat2str(seq{k}.desc.amp) ',''z'',1,NOWAIT );\n']);
				 fprintf(saida,['delay(' mat2str(seq{k}.desc.temp) ');\n']);
				  
				  else
				  fprintf(saida,['zgradpulse(' mat2str(seq{k}.desc.amp) ',' ...
					  mat2str(seq{k}.desc.temp) ');\n']);
				  
				  end
				  fprintf(saida,['delay(0.01);\n']);
   
    end
      
end

 fprintf(saida,'\n');
 fprintf(saida,'} \n'); 
 fclose all;
  

 
 unix(['seqgen ' nome]);


 
 
 
 
 
 
 
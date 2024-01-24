function pul = writepul(pul,seq,seq2)

locca = pwd; 
cd ..;
cd ..;
cd shapes;
n = length(pul)-1;
d =1;
ddir2 = '/home/quantum/vnmrsys/shapelib/'; 
%ddir2 = [];
for p=1:n
    switch pul{p}.tip
        
        case 'pulse'
              
              clear ang phase targ shp;
             for k=1:length(pul{p}.desc.forma)
                  
                  chan = pul{p}.desc.chan{k};
                  
                  for m =1:length(pul{p}.desc.forma{k})
                      
                  shp{m} = load([pul{p}.desc.forma{k}{m}]);
                  
                  ang(m) = pul{p}.desc.ang(m);
                  
                  phase(m) = pul{p}.phase(m);
                  
                  targ(m) = pul{p}.desc.targ(m);
              
                  end
                                    
                  typs = pul{p}.desc.typs{k};
                  
                  [wr{k} ph{k} pul{p}.desc.power(k) np dt] = ...
                  composepul(shp, pul{p}.desc.pw,ang,phase,targ,typs);
                  
                  r  = 1023*wr{k}/max(wr{k});
              
                  th = mod(ph{k}*180/pi,360); 
                  
                  pul{p}.desc.formaspec{k} = [seq2 'pul' mat2str(d) 'Ch' mat2str(chan) '.RF'];
                  
                  saida = fopen( pul{p}.desc.formaspec{k} ,'w'); 
                   
                  for n = 1:length(r)
                   fprintf(saida,'%10.6f %12.3f %12.6f\n', th(n),r(n),1.0);
                  end
                  

                 fclose all;
                if isunix
                 unix(['cp ' pul{p}.desc.formaspec{k} '  ' ddir2 pul{p}.desc.formaspec{k}]);
                end
             end
    d =d+1;
    end
  
end

eval(['cd ' locca])
    end

    
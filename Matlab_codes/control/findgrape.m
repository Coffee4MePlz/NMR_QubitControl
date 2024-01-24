function  findgrape(mol,spectro,grape)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initial set up
grape.nch = length(grape.chanel);  

for m=1:grape.nch; grape.h{2*m-1} = mol.Ich{grape.nch(m),1}; grape.h{2*m} = mol.Ich{grape.nch(m),2}; end

grape.H0{1}  =  mol.Hzee  +  mol.Hint; n =1;

for k=1:mol.nspin;  
    for m= 1:length(grape.chem{k});  
        n= n+1; grape.H0{n} = grape.H0{1} + 2*pi* grape.chem{k}(m) * mol.Iz{k}; 
    end
end

grape.ndq = length(grape.H0);

for k=1:grape.ndq
            ur1 = expm(-1i*grape.H0{k}*spectro.rof1*1e-6);
            ur2 = expm(-1i*grape.H0{k}*spectro.rof2*1e-6);
            grape.Utarg{k} = (ur2')*grape.U*(ur1');
end


grape.N = length(grape.H0{1});

grape.nrf = length(grape.rf);
 
if grape.typ == 'c'

    options=optimset('Diagnostics','on','Gradobj','on','Hessian','off',...
 'HessUpdate','bfgs','MaxFunEvals',1e10,'TolFun',1e-8,'TolX',1e-14,'TolCon',1e-14,...
 'MaxIter',1500,'Display','Iter','LargeScale','off','Algorithm','interior-point');
    

else
    
    options=optimset('Diagnostics','on','Gradobj','on','Hessian','off',...
 'HessUpdate','bfgs','MaxFunEvals',1e10,'TolFun',1e-14,'TolX',1e-14,...
 'MaxIter',1500,'Display','Iter','LargeScale','off');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% with optimization with 20% of the points

x = [];

for k=0.2:0.8:1.0
    
     np = round(grape.np*k);
    
     [x grape]=smoothg(grape,np,x);
    
    if grape.typ == 'c'
    
      y = linspace(1,np,np)/np;
      y = 2*(y - mean(y));  
      y = -cosh(y*10);  
      y = (y - min(y));  
      y =2*pi*grape.maxpower*y/max(y);
      penal= repmat(y,1,grape.nch*2)';
      [x, f] = fmincon(@(x) ffun(x,grape),x,[],[],[],[],-penal,penal,[],options); 
     
    else
      [x f] = fminunc(@(x) ffun(x,grape),x,options);
    end
    
end


disp(['SOLUTION WITH FIDELITY ' mat2str(1-f)]);


%-----------------------------------------------------------------------
% saving solution
for k=1:grape.nch
       
       n = 1 + (2*grape.np2*(k-1));
       
       hx = x(n:(n+grape.np2-1));
       
       hy = x((n + grape.np2):(n+2*grape.np2-1));
       
       figure(k)
      
       plot(hx,'-r'); hold on; plot(hy,'-k');
 
       [th{k}, rr{k}] = cart2pol(hx,hy);
 
 
       th{k} = mod(180*th{k}/pi,360);
 

       if  grape.chanel(k) == 1;  
              saida = fopen([grape.nome 'Obs.RF'],'w'); 
       else
              saida = fopen([grape.nome 'Dec.RF'],'w'); 
       end 
 
       
    for m = 1:length(rr{k})
        fprintf(saida,'%10.6f %12.3f %12.6f\n', th{k}(m),rr{k}(m),1.0);
    end
 
    fclose all;
 end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [f, g] = ffun(x1,grap)
    %function for fidelity and gradient        
    
    f = 0; g = 0;  
   
    aux2 = zeros(grap.np2*2*grap.nch,1);
    Uf = cell(1,grap.np2);
    Ub = cell(1,grap.np2);
    
    
    for p = 1:grap.ndq;
    
        for r =1:grap.nrf;
           
            xrf = x1*(1+grap.rf(r)) ;        
            
            U = eye(grap.N);
            
            for k=1:grap.np2;
           
                H = grap.H0{p};
                
                for m=1:(2*grap.nch); H = H + xrf(k + grap.np2*(m-1)) * grap.h{m}; end
                
                Uf{k} = expm(-1i*grap.dt*H)*U; U = Uf{k}; %forward;
                
            end
            
            for k=1:grap.np2; Ub{k} = Uf{k}*Uf{end}'; end %backward
             
            
            %--------fidelity ------------------
            
            aux = 1 - trace(grap.Utarg{p}' * U) * trace(U'*grap.Utarg{p})/ (grap.N)^2 ;
             
            f = f + aux;
            
            %---------- gradient --------------
          
            for k=1:grap.np2
                    
                   Pk = Ub{k}*grap.Utarg{p};
                   Xk = Uf{k};
                   A = -1i*grap.dt * trace(Xk'*Pk) ;
                   
                   for m=1:(2*grap.nch)
              
                       
                       aux2(k + grap.np2*(m-1)) =  real(trace(Pk' * (grap.h{m}) * Xk) * A);
                    
                   end
               
            end
          
            g = g  + aux2 ;         
        
        end
        
    end
    f = f/(grap.nrf*grap.ndq);
    g = -2*g/(grap.nrf*grap.ndq*grap.N^2); 
    
    end
    
function [x0 grap]=smoothg(grap,nf,x);
    
  if isempty(x)
        x0 = [];
        for m=1:grap.nch
            x0 = [x0 2*pi*grap.maxpower *( 0.1*(rand(1,nf)-0.5))];
            x0 = [x0 2*pi*grap.maxpower *( 0.1*(rand(1,nf)-0.5))];
        end
        
        grap.dt = grap.temp/nf;
    else
        
     ni = length(x)/(2*grap.nch);   
     xt2 = [0 (1:ni)*(grap.dt)-grap.dt/2 grap.temp];
     grap.dt = grap.temp/nf;
     xt = (0:nf)*(grap.dt);
     x0 = [];
      for m=1:(2*grap.nch)    
              n =  grap.np2*(m-1);
              h = x((n+1):(n + grap.np2));
              x0 = [x0 (spline(xt2,[0 h' 0],xt))];
      end
    
  end
 grap.np2 = length(x0)/(2*grap.nch); ;
 x0 = x0';
end
        
   
 

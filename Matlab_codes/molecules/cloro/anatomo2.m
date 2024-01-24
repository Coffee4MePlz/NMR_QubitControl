clear all
spectropar;
cloroformio;
[tran1 tran2 ff1 opr] = transitions(mol,1);


spechr= varian2mat('refh',64*1024);
spechr= apodize(spechr);
spechr{1} = fitespec(spechr{1},ff1,0.01,[0.6 0.1],[-200 200]);
plotfit(1,spechr{1},[-200 200]);


speccr= varian2mat('refc',64*1024);
speccr= apodize(speccr);
speccr{1} = fitespec(speccr{1},ff1,0.58,[0.6 0.1],[-200 200]);
plotfit(1,speccr{1},[-200 200]);


nh = abs(spechr{1}.ILP(1) + spechr{1}.ILP(3))/2;
nc = abs(speccr{1}.ILP(1) + speccr{1}.ILP(3))/2;
f =4*nc/nh;

%f = 0.7355;

aa=varian2mat('J0Hm',64*1024);
bb=varian2mat('J0Cm',64*1024);
aa = apodize(aa);
bb = apodize(bb);

for k=1:length(aa);
	aa{k} =  fitespec(aa{k},ff1,0.01,[0.6 0.1],[-200 200]);
    plotfit(1,aa{k},[-200 200]); pause(0.2);
	
	bb{k} =  fitespec(bb{k},ff1,0.58,[0.6 0.1],[-200 200]);
    plotfit(1,bb{k},[-200 200]); pause(0.2);
end

for m=1:61
for k=1:4;
    struc{m}{1,k} =  f*aa{4*(m-1) + k}.ILP;
    struc{m}{2,k} =  bb{4*(m-1) + k}.ILP;
end

mx(m) = ( bb{4*(m-1) + 1}.ILP(1) +  bb{4*(m-1) + 1}.ILP(3) )/(7.7171e5);
my(m) = ( bb{4*(m-1) + 1}.ILP(2) +  bb{4*(m-1) + 1}.ILP(4) )/(7.7171e5);
end
close all

t =  (0:60)*4e-4;
plot(t,mx,'or',t,my,'ob');

for m=1:61
roh{m} = tomography({'+II','+IX','+IY','+XX'},struc{m},2);
roh{m} = eye(4)/4 + (0.25/3.8794e+05)*roh{m};
showmat(2,roh{m}); pause(0.1);
end

save J0menos  roh t mx my

% 
% clc
% 
% roh{m} = eye(4)/4 + (0.25/3.8794e+05)*roh{m};
% showmat(1,roh);
% %roh2 = diag([0.75 -0.25 -0.25 -0.25]);
% roh2 = diag([1 0 0 0]); 
% rr = cos(pi/3)*[1;0] + sin(pi/3)*[0;1]; rr = diag(diag(rr));
% roh2 = kron(rr*rr',[1 1;1 1]/2);
% 
% showmat(2,roh2);
% fidelmat(roh2,roh)
% 
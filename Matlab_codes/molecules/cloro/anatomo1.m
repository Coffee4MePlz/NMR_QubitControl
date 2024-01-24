clear all
spectropar;
cloroformio;
[tran1 tran2 ff1 opr] = transitions(mol,1);


spechr= varian2mat('refh',16*1024);
spechr{1} = fitespec(spechr{1},ff1,0.58,[0.6 3],[-200 200]);

speccr= varian2mat('refc',16*1024);
speccr{1} = fitespec(speccr{1},ff1,0.58,[0.6 3],[-200 200]);

nh = abs(spechr{1}.ILP(1) + spechr{1}.ILP(3))/2;
nc = abs(speccr{1}.ILP(1) + speccr{1}.ILP(3))/2;
f =4*nc/nh;

f = 0.0294;
aa=varian2mat('bd3h',64*1024);
bb=varian2mat('bd3c',64*1024);
aa = apodize(aa);
bb = apodize(bb);

for k=1:length(aa);
	aa{k} =  fitespec(aa{k},ff1,0.58,[0.6 3],[-200 200]);
    plotfit(1,aa{k},[-200 200]); pause(0.2);
	
	bb{k} =  fitespec(bb{k},ff1,0.58,[0.6 3],[-200 200]);
    plotfit(1,bb{k},[-200 200]); pause(0.2);
end

for m=1:51;
	for k=1:4;
    struc{m}{1,k} =  f*aa{4*(m-1) + k}.ILP;
    struc{m}{2,k} =  bb{4*(m-1) + k}.ILP;
	end
	mx(m) = bb{4*(m-1) + 1}.ILP(1) + bb{4*(m-1) + 1}.ILP(3);
	mxb(m) = f*aa{4*(m-1) + 1}.ILP(1) + f*aa{4*(m-1) + 1}.ILP(3);
end
mx2 = mx/6.5742e+04;
mx1 = mxb/6.5742e+04;

for m=1:51; 
	roh{m} = tomography({'+II','+IX','+IY','+XX'},struc{m},2);
	roh{m} = eye(4)/4 + (0.25/2.4280e+04)*roh{m};
end
t  =(0:51)*(2/215);

% t = linspace(0,0.02,51);
% rr = cos(pi/3)*[1;0] + sin(pi/3)*[0;1];  
% roh2 = kron(rr*rr',[1 -1;-1 1]/2);
% for k=1:51;
% 	 U = zzgate(2,2*pi*mol.j(1,2)*t(k),1,2);
% 	 roh3 = U*roh2*U' ;
%      ff(k) = fidelmat(roh{k},roh3);
% 	 xx(k) = trace(roh3*multigate(2,2,xgate));
% end
% plot(t,ff) ;figure; plot(t,xx,'ob',t,mx2,'sr')

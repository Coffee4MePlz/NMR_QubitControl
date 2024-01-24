clear all
spectropar;
cloroformio;
[tran1 tran2 ff1 opr] = transitions(mol,1);


spechr= varian2mat('refhn',64*1024);
spechr= apodize(spechr);
spechr{1} = fitespec(spechr{1},ff1,0.01,[0.6 0.1],[-200 200]);
plotfit(1,spechr{1},[-200 200]);


speccr= varian2mat('refcn',64*1024);
speccr= apodize(speccr);
speccr{1} = fitespec(speccr{1},ff1,0.58,[0.6 0.1],[-200 200]);
plotfit(1,speccr{1},[-200 200]);


nh = abs(spechr{1}.ILP(1) + spechr{1}.ILP(3))/2;
nc = abs(speccr{1}.ILP(1) + speccr{1}.ILP(3))/2;
f =4*nc/nh;
 f = 0.0294;
%f = 0.7355;

aa=varian2mat('ppsh3',64*1024);
bb=varian2mat('ppsc3',64*1024);
aa = apodize(aa);
bb = apodize(bb);

for k=1:length(aa);
	aa{k} =  fitespec(aa{k},ff1,0.58,[0.6 0.1],[-200 200]);
    plotfit(1,aa{k},[-200 200]); pause(0.2);
	
	bb{k} =  fitespec(bb{k},ff1,0.58,[0.6 0.1],[-200 200]);
    plotfit(1,bb{k},[-200 200]); pause(0.2);
end


for k=1:4;
    struc{1,k} =  f*aa{k}.ILP;
    struc{2,k} =  bb{k}.ILP;
end


roh = tomography({'+II','+IX','+IY','+XX'},struc,2);
clc

%roh = eye(4)/4 + (0.25/2.4280e+04)*roh;
showmat(1,roh);
roh2 = diag([0.75 -0.25 -0.25 -0.25]);
roh2 = diag([1 0 0 0]); 

rr1=[0.425 -0.147+0.101i -0.147+0.101i 0.425;
-0.147-0.101i 0.075 0.075 -0.147-0.101i;
-0.147-0.101i 0.075 0.075 -0.147-0.101i;
0.425 -0.147+0.101i -0.147+0.101i 0.425];


rr2 = [0.075 -0.179 -0.179 0.075;
-0.179 0.425 0.425 -0.179;
-0.179 0.425 0.425 -0.179;
0.075 -0.179 -0.179 0.075];

rr3=[-0.437 -0.368+0.252i -0.368+0.252i 0.187;
-0.368-0.252i 0.437 1.006 -0.368-0.252i ;
-0.368-0.252i 1.06 0.437  -0.368-0.252i ;
0.187 -0.368+0.252i -0.368+0.252i -0.437];


%rr = cos(pi/3)*[1;0] + sin(pi/3)*[0;1]; rr = diag(diag(rr));
%roh2 = kron(rr*rr',[1 1;1 1]/2);

showmat(2,rr3);
fidelmat(roh,rr3)


clc
clear all
trifluor1;
[freq1 tran ff1 opR opI] = transitions(mol,1);
ff1 = ff1+0.5;

trifluor2;
[freq2 tran ff2 opR opI] = transitions(mol,2);
ff2 = ff2+0.5;


trifluor3;

[freq3 tran ff3 opR opI] = transitions(mol,3);
ff3 = ff3+0.5;


sss1 = {'est000n1','est000n2','est000n3','est000n4','est000n5','est000n6','est000n7'};
sss2 = {'est0002n1','est0002n2','est0002n3','est0002n4','est0002n5','est0002n6','est0002n7'};
sss3 = {'est0003n1','est0003n2','est0003n3','est0003n4','est0003n5','est0003n6','est0003n7'};



for k=1:7
spec1{k}= varian2mat(sss1{k},32*1024);

spec2{k}= varian2mat(sss2{k},32*1024);

spec3{k}= varian2mat(sss3{k},32*1024);

spec1{k} = apodize(spec1{k}); spec2{k} = apodize(spec2{k});spec3{k} = apodize(spec3{k});
end

for k=1:7


	spec1{k}{1} = fitespec(spec1{k}{1},ff1,T2(1),[0.1 2],[-400 400]);
    
    spec2{k}{1} = fitespec(spec2{k}{1},ff2,mol.T2(2),[0.1 2],[-400 400]);

	spec3{k}{1} = fitespec(spec3{k}{1},ff3,T2(3),[0.1 2],[-400 400]);
     
	plotfit(k,spec1{k}{1},[-100 100]); pause(0.2)
    plotfit(k,spec2{k}{1},[-200 200]); pause(0.2)
	plotfit(k+10,spec3{k}{1},[-150 150]); pause(0.2)
end

for k=1:7;
    struc{1,k} =  spec1{k}{1}.ILP;
    struc{2,k} =  spec1{k}{1}.ILP;
    struc{3,k} =  spec3{k}{1}.ILP;
end

%--------------------------------------------------------------------------

roh = tomography({'+III','+IIY','+IYY','+YII','+XXY','+XYX','+XXX'},struc,3);
% roh = 0.1925*roh /1.0976e+06 + eye(4)/4
showmat(30,roh);

% rr= roty(140*pi/180)*diag([1 0])*roty(-140*pi/180); rr= diag(diag(rr));
% roh2 = kron(eye(2)/2,rr) ;
%  ux = multigate(2,[2],rotx(60*pi/180));
% % uxb = multigate(2,[1 2],rotx(-pi/2));
%  uy = multigate(2,[2],roty(-60*pi/180));
% % uyb = multigate(2,[1 2],roty(-pi/2));
%  zz = zzgate(2,pi,1,2);
%  u = uy*zz*ux;
% %u = controlgate(2,1,2,xgate);
% roh2 = u*roh2*u';
% showmat(40,roh2)
% fidelmat(roh,roh2)





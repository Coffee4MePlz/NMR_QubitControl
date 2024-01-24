clc

trifluor;

[freq1 tran ff1 opR opI] = transitions(mol,2);


%sss1 = {'Iii','Ixi','Iyi','Iix','Iiy','Ixx','Iyy','Ixy','Iyx'};
%sss1 = {'0ii','0xi','0yi','0ix','0iy','0xx','0yy','0xy','0yx'};
%sss1 = {'30ii','30xi','30yi','30ix','30iy','30xx','30yy','30xy','30yx'};
sss1 = {'1802ii','1802xi','1802yi','1802ix','1802iy','1802xx','1802yy','1802xy','1802yx'};
%sss1 = {'90ii','90xi','90yi','90ix','90iy','90xx','90yy','90xy','90yx'};
%sss1 = {'120ii','120xi','120yi','120ix','120iy','120xx','120yy','120xy','120yx'};
%sss1 = {'150ii','150xi','150yi','150ix','150iy','150xx','150yy','150xy','150yx'};
%sss1 = {'180ii','180xi','180yi','180ix','180iy','180xx','180yy','180xy','180yx'};

for k=1:9
spec1{k}= varian2mat(sss1{k},32*1024);
spec1{k} = apodize(spec1{k}); 
end

for k=1:9
spec1{k}{1} = fitespec(spec1{k}{1},ff1,mol.T2(2),[0.01 1],[-400 400]);
plotfit(k,spec1{k}{1},[-200 200]); pause(0.2)
end
bb = [1 4 5];
for k=1:9;
%	if any(k==bb); 
%		aa = spec1{k}{1}.ILP; 
%		spec1{k}{1}.ILP([1 2 3 4]) = aa([5 6 7 8]); 
%		spec1{k}{1}.ILP([5 6 7 8]) = aa([1 2 3 4]); 
%	end
    struc{1,k} =  spec1{k}{1}.ILP;
end


roh = tomography2({'+III','+XII','+YII','+IIX','+IIY','+XIX','+YIY','+XIY','+YIX'},...
struc,3,[1 3],2);

%roh = tomography2({'+III','+IXI','+IYI','+IIX','+IIY','+IXX','+IYY','+IXY','+IYX'},...
%struc,3,[2 3],1);

clc

roh = roh/6.102796299574417e+05;

disp(roh)

showmat(30,roh);

roh2 = kron(eye(2)/2,diag([1 0])) ;
u = multigate(2,2,roty(40*pi/180)); roh2 = diag(diag(u*roh2*u'));
%u = controlgate(2,1,2,xgate); 
%u = swapgate(2,1,2);  roh2 = u*roh2*u';
u = multigate(2,[1 2],roty(pi/2)); u = zzgate(2,180*pi/180,1,2)*u;
%u = multigate(2,[1 2],rotz(-pi/2))*u;
u = multigate(2,[1 2],roty(-pi/2))*u; u = zzgate(2,180*pi/180,1,2)*u;
u = multigate(2,[1 2],rotx(pi/2))*u; u = zzgate(2,180*pi/180,1,2)*u;
u = multigate(2,[1 2],rotx(-pi/2))*u;
%u = multigate(2,[1 2],rotz(pi))*u;

roh2 = u*roh2*u';

disp(roh2)
showmat(31,roh2);
fidelmat(roh,roh2)

trace(roh)
trace((roh)^2)
trace(roh2^2)





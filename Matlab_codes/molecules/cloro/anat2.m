clear all
spectropar;
cloroformio;
[tran1 tran2 ff1 opr] = transitions(mol,1);


t = linspace(0,420,15);
%t = linspace(0,10,20);
%t = linspace(0,20,20);
%t = linspace(0,50,20);
aa=varian2mat('t2',8*1024);


aa = apodize(aa);

for k=1:length(aa);
	aa{k} =  fitespec(aa{k},ff1-5,0.58,[0.6 3],[-200 200]);
    plotfit(1,aa{k},[-200 200]); pause(0.2);
	mx(k) = sum(aa{k}.ILP(1:2:3));
	my(k) = sum(aa{k}.ILP(2:2:4)); 
	mm(k) = sqrt(mx(k)^2 + my(k)^2);
end

figure(10)
plot(t,mm,'ok') ; 
hold on
fun{1} = @(c,x) exp(-x/c(1));

[INLP ILP] = fminspleas(fun,[200],t,mm,0,400);

xx = linspace(0,max(t),100);
yy = ILP*fun{1}(INLP,xx);
plot(xx,yy,'r')
disp(INLP)

figure(11)
plot(t,mx,'ok',t,my,'or') ;
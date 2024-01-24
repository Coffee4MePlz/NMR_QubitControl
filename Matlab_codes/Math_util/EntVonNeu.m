function Evn = EntVonNeu(roh);

autoval = real(eig(roh));

S = 0;

for k = 1:length(autoval);
    if autoval(k,1) > 0
       S = -(autoval(k,1))*log2(autoval(k,1)) + S;
    end
end

Evn = S;
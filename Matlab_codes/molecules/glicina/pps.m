

 spectropar;

 glicina;
 
 roh = 4*(mol.Iz{1} + mol.Iz{2}) + mol.Iz{3}+ mol.Iz{4};
 
 pw = 7;
 
 [U roh] = decshapedpulse('isech200',50000,0.0,-170.4,10,roh)
 
 [U roh] = decrgpulse(pw,90,90,roh);
 
 [spec]  = specgen(roh,2);
 
 plotfid(1,spec);
 
 plotesp(2,spec,[-240 240])
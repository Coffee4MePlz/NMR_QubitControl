

 spectropar;

 tricloro;
 
 roh = 4*mol.Iz{1} + mol.Iz{2} + mol.Iz{3};
 
 pw = 7;

 [U roh] = decrgpulse(pw,90,90,roh);
 
 [spec]  = specgen(roh,2);
 
 plotfid(1,spec);
 
 plotesp(2,spec,[-800 800])
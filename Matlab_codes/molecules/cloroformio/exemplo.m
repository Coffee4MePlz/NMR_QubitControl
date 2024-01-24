

 spectropar;

 cloroformio;
 
 roh = 4*mol.Iz{1} + mol.Iz{2};
 
 pw = 7;
 
 p1 = 10;
 
 d2 =abs( 1/(2*mol.j(1,2)));
 
 [U roh] = rgpulse(pw,90,0,roh);
 
 [U roh] =delay(d2,roh);
 
 %[U roh] = decrgpulse(p1/2,45,0,roh);
 
 [spec]  = specgen(roh,1);
 
 %plotfid(1,spec);
 
 plotesp(5,spec,[-200 200])
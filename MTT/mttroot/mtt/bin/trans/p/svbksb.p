PROCEDURE svbksb(u: glmpbynp; w: glnparray; v: glnpbynp;
       m,n: integer; b: glmparray; VAR x: glnparray);
{

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Version control history
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % $Id: svbksb.p,v 1.3 1998/08/12 11:28:13 peterg Exp $
% % $Log: svbksb.p,v $
% % Revision 1.3  1998/08/12 11:28:13  peterg
% % Zapped unnecessary np and mp arguments
% %
% % Revision 1.2  1998/08/12 11:09:54  peterg
% % Taken from NR share library nrpas13 as  SVBKSB.PAS
% % and renamed  svbksb.p
% %
% % Revision 1.1  1998/08/12 11:09:02  peterg
% % Initial revision
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 }


(* Programs using SVBKSB must define the types
TYPE
   glnparray = ARRAY [1..np] OF real;
   glmparray = ARRAY [1..mp] OF real;
   glnpbynp = ARRAY [1..np,1..np] OF real;
   glmpbynp = ARRAY [1..mp,1..np] OF real;
in the main routine. *)
VAR
   jj,j,i: integer;
   s: real;
   tmp: glnparray;
BEGIN
   FOR j := 1 TO n DO BEGIN
      s := 0.0;
      IF (w[j] <> 0.0) THEN BEGIN
         FOR i := 1 TO m DO BEGIN
            s := s+u[i,j]*b[i]
         END;
         s := s/w[j]
      END;
      tmp[j] := s
   END;
   FOR j := 1 TO n DO BEGIN
      s := 0.0;
      FOR jj := 1 TO n DO BEGIN
         s := s+v[j,jj]*tmp[jj];
      END;
      x[j] := s
   END
END;

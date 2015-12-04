PROCEDURE mtt_ludcmp(VAR a: StateMatrix; n: integer;
		     VAR indx: StateVector; VAR d: real);

{ This file is derived from the NUMERICAL RECIPES PASCAL SHAREWARE DISKETTE.
 Please read the README file in $MTTPATH/trans/p
 }

{
###############################################################
## Version control history
###############################################################
## $Id: mtt_ludcmp.p,v 1.2 1998/08/17 15:52:54 peterg Exp $
## $Log: mtt_ludcmp.p,v $
## Revision 1.2  1998/08/17 15:52:54  peterg
## Added notice.
##
## Revision 1.1  1998/08/14 07:58:13  peterg
## Initial revision
##
###############################################################
}

(* Programs using LUDCMP must define the types
TYPE
   glnpbynp = ARRAY [1..np,1..np] OF real;
   glnarray = ARRAY [1..n] OF real;
   glindx = ARRAY [1..n] OF integer;
in the main routine. *)

CONST
   tiny=1.0e-20;

VAR
   k,j,imax,i: integer;
   sum,dum,big: real;
   vv: glnarray;

BEGIN
   d := 1.0;
   FOR i := 1 TO n DO BEGIN
      big := 0.0;
      FOR j := 1 TO n DO IF (abs(a[i,j]) > big) THEN big := abs(a[i,j]);
      IF (big = 0.0) THEN BEGIN
         writeln('pause in LUDCMP - singular matrix'); readln
      END;
      vv[i] := 1.0/big
   END;
   FOR j := 1 TO n DO BEGIN
      FOR i := 1 TO j-1 DO BEGIN
         sum := a[i,j];
         FOR k := 1 TO i-1 DO BEGIN
            sum := sum-a[i,k]*a[k,j]
         END;
         a[i,j] := sum
      END;
      big := 0.0;
      FOR i := j TO n DO BEGIN
         sum := a[i,j];
         FOR k := 1 TO j-1 DO BEGIN
            sum := sum-a[i,k]*a[k,j]
         END;
         a[i,j] := sum;
         dum := vv[i]*abs(sum);
         IF (dum > big) THEN BEGIN
            big := dum;
            imax := i
         END
      END;
      IF (j <> imax) THEN BEGIN
         FOR k := 1 TO n DO BEGIN
            dum := a[imax,k];
            a[imax,k] := a[j,k];
            a[j,k] := dum
         END;
         d := -d;
         vv[imax] := vv[j]
      END;
      indx[j] := imax;
      IF (a[j,j] = 0.0) THEN a[j,j] := tiny;
      IF (j <> n) THEN BEGIN
         dum := 1.0/a[j,j];
         FOR i := j+1 TO n DO BEGIN
            a[i,j] := a[i,j]*dum
         END
      END
   END;
END;
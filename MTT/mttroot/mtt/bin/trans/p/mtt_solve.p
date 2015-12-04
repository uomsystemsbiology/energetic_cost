PROCEDURE mtt_solve(VAR x : StateVector;
			A : StateMatrix;
			B : StateVector;
			n : integer);

{

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Version control history
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % $Id: mtt_solve.p,v 1.4 2000/08/01 09:33:28 peterg Exp $
% % $Log: mtt_solve.p,v $
% % Revision 1.4  2000/08/01 09:33:28  peterg
% % Include mtt_solve_lud within scope of this function - seems to avaoid
% % probs with p2c
% %
% % Revision 1.3  1999/10/26 23:44:36  peterg
% % Now packages up mtt_solve_lud OR mtt_solve_svd
% %
% % Revision 1.2  1998/08/14 12:09:13  peterg
% % A passed by value - its destroyed by SVDcm
% %
% % Revision 1.1  1998/08/13 08:51:57  peterg
% % Initial revision
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 }

VAR
   i	     : integer;
   wmax,wmin : real;
   w	     : StateVector ;
   v	     : StateMatrix;


(*$I $MTTPATH/trans/p/mtt_solve_lud.p *)
(* I $MTTPATH/trans/p/mtt_solve_svd.p *)

BEGIN{mtt_solve}
      mtt_solve_lud(x,A,B,N);
      (*** mtt_solve_svd(x,A,B,N); ***)
END{mtt_solve};

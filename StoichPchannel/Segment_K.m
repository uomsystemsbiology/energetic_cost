function [K i_lin] = Segment_K(mttpar)
## function [K i_lin] = Segment_K(mttpar)
## Generated by odeo2kpar on Mon Mar 30 09:29:40 AEDT 2015

## Parameters 
   c_m 	= mttpar(1);
   f 	= mttpar(2);
   k_c 	= mttpar(3);
   k_ch 	= mttpar(4);
   k_cm 	= mttpar(5);
   k_o 	= mttpar(6);
   k_oh 	= mttpar(7);
   k_om 	= mttpar(8);
   k_k 	= mttpar(9);
   k_l 	= mttpar(10);
   k_n 	= mttpar(11);
   rt 	= mttpar(12);
   v 	= mttpar(13);
   v_eq 	= mttpar(14);
   g_s 	= mttpar(15);
   kappa_k 	= mttpar(16);
   kappa_l 	= mttpar(17);
   kappa_n 	= mttpar(18);
   v_k 	= mttpar(19);
   v_l 	= mttpar(20);
   v_n 	= mttpar(21);
   v_unit 	= mttpar(22);
   x_g 	= mttpar(23);
   z 	= mttpar(24);

## Free-energy constants
  K(1,1) =  f**2/(c_m*rt);
  K(2,1) =  k_c;
  K(3,1) =  k_o;
  K(4,1) =  k_k;
  K(5,1) =  k_k;
  K(6,1) =  k_cm;
  K(7,1) =  k_om;
  K(8,1) =  k_ch;
  K(9,1) =  k_oh;
  K(10,1) =  k_n;
  K(11,1) =  k_n;
  K(12,1) =  k_l;
  K(13,1) =  k_l;
if (min(K)==0); warning("K has a zero element"); endif
## Index of linear Cs
i_lin = [];
i_lin = [i_lin 1 ];
endfunction

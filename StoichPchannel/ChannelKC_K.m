function [K i_lin] = ChannelKC_K(mttpar)
## function [K i_lin] = ChannelKC_K(mttpar)
## Generated by odeo2kpar on Thu Jul  9 10:27:21 AEST 2015

## Parameters 
   k_c 	= mttpar(1);
   k_o 	= mttpar(2);
   k_k 	= mttpar(3);
   rt 	= mttpar(4);
   v_e 	= mttpar(5);
   v_i 	= mttpar(6);
   kappa_k 	= mttpar(7);
   v_k 	= mttpar(8);
   x_g 	= mttpar(9);
   z 	= mttpar(10);

## Free-energy constants
  K(2,1) =  k_c;
  K(3,1) =  k_o;
  K(4,1) =  k_k/v_i;
  K(5,1) =  k_k/v_e;
if (min(K)==0); warning("K has a zero element"); endif
## Index of linear Cs
i_lin = [];
endfunction
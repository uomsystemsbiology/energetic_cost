
function [mtty] = ChannelNaC_odeo(mttx,mttu,mttt,mttpar);
## [mtty] = ChannelNaC_odeo(mttx,mttu,mttt,mttpar);
## System ChannelNaC, representation odeo, language m; 
## File ChannelNaC_odeo.m; 
## Generated by MTT on Fri Dec  4 08:42:27 AEDT 2015; 

## Horrible fudge to make mtt_m2p work
global ...
mtt_no_globals ;

## Parameters 
   k_ch 	= mttpar(1);
   k_cm 	= mttpar(2);
   k_oh 	= mttpar(3);
   k_om 	= mttpar(4);
   k_n 	= mttpar(5);
   rt 	= mttpar(6);
   v_e 	= mttpar(7);
   v_i 	= mttpar(8);
   kappa_n 	= mttpar(9);
   v_n 	= mttpar(10);
   x_g 	= mttpar(11);
   z 	= mttpar(12);

## States 
mttx1 = mttx(1) ;
mttx2 = mttx(2) ;
mttx3 = mttx(3) ;
mttx4 = mttx(4) ;
mttx5 = mttx(5) ;
mttx6 = mttx(6) ;
mttx7 = mttx(7) ;

## Inputs 
mttu1 = mttu(1);
mttu2 = mttu(2);
mttu3 = mttu(3);

## Unknown Inputs 

   mtty = zeros(3,1);

## BEGIN Code

## Code
  mtty(1) = mma_cr(2,mttx1)*(e**((3*mttx1)/rt)*k_cm*mttx2 - k_om*mttx3);
  mtty(2) = mma_cr(3,mttx1)*(e**((4*mttx1)/rt)*k_ch*mttx4 - k_oh*mttx5);
  mtty(3) = (ghk_fun(mttx1,v_n)*k_ch*k_n*k_om**3*kappa_n*mttx3**3*mttx4*(e**(( \
  mttx1*z)/rt)*mttx6*v_e - mttx7*v_i))/(v_e*v_i*x_g**4);
## END Code
endfunction

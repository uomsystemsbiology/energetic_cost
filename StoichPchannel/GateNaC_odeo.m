
function [mtty] = GateNaC_odeo(mttx,mttu,mttt,mttpar);
## [mtty] = GateNaC_odeo(mttx,mttu,mttt,mttpar);
## System GateNaC, representation odeo, language m; 
## File GateNaC_odeo.m; 
## Generated by MTT on Wed Mar 11 09:48:23 AEDT 2015; 

## Horrible fudge to make mtt_m2p work
global ...
mtt_no_globals ;

## Parameters 
   k_ch 	= mttpar(1);
   k_cm 	= mttpar(2);
   k_oh 	= mttpar(3);
   k_om 	= mttpar(4);
   rt 	= mttpar(5);

## States 
mttx1 = mttx(1) ;
mttx2 = mttx(2) ;
mttx3 = mttx(3) ;
mttx4 = mttx(4) ;
mttx5 = mttx(5) ;

## Inputs 
mttu1 = mttu(1);
mttu2 = mttu(2);

## Unknown Inputs 

   mtty = zeros(3,1);

## BEGIN Code

## Code
  mtty(1) = rt*(log(k_ch*mttx4) + 3*log(k_om*mttx3) - log(k_ch) - 3*log(k_om));
  mtty(2) = mma_cr(2,mttx1)*(e**((3*mttx1)/rt)*k_cm*mttx2 - k_om*mttx3);
  mtty(3) = mma_cr(3,mttx1)*(e**((4*mttx1)/rt)*k_ch*mttx4 - k_oh*mttx5);
## END Code
endfunction

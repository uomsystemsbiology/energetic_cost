
function [mttpar] = ChannelNaC_numpar;
## [mttpar] = ChannelNaC_numpar;
## System ChannelNaC, representation numpar, language m; 
## File ChannelNaC_numpar.m; 
## Generated by MTT on Fri Dec  4 08:42:27 AEDT 2015; 

## Horrible fudge to make mtt_m2p work
global ...
mtt_no_globals ;

   mttpar = zeros(12,1);

## BEGIN Code

## User defined code from ChannelNaC_numpar.txt
  v_i =  1;
  v_e =  1;
  k_n =  1/1000;
  rt =  8.3144621*300;		
  kappa_n =   0.13204e-9;		
  z	=  1;			
  v_n	=  0;
  k_cm =   105.49;
  k_om =   1;
  k_ch =   1;
  k_oh =     6.3281e-05;
  x_g  =  0.001;

  ## Set up the mttpar vector
  mttpar(1)	= k_ch;
  mttpar(2)	= k_cm;
  mttpar(3)	= k_oh;
  mttpar(4)	= k_om;
  mttpar(5)	= k_n;
  mttpar(6)	= rt;
  mttpar(7)	= v_e;
  mttpar(8)	= v_i;
  mttpar(9)	= kappa_n;
  mttpar(10)	= v_n;
  mttpar(11)	= x_g;
  mttpar(12)	= z;
## END Code
endfunction

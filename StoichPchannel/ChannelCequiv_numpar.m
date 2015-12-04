
function [mttpar] = ChannelCequiv_numpar;
## [mttpar] = ChannelCequiv_numpar;
## System ChannelCequiv, representation numpar, language m; 
## File ChannelCequiv_numpar.m; 
## Generated by MTT on Mon Mar  9 08:37:01 AEDT 2015; 

## Horrible fudge to make mtt_m2p work
global ...
mtt_no_globals ;

   mttpar = zeros(6,1);

## BEGIN Code

## User defined code from ChannelCequiv_numpar.txt
  k	=  1.0; 
  rt	=  8.3144621*300;
  g	=  1.0;
  kappa	=  1.0;
  z	=  2;
  k_e	=  1; 

  ## Set up the mttpar vector
  mttpar(1)	= k;
  mttpar(2)	= k_e;
  mttpar(3)	= rt;
  mttpar(4)	= g;
  mttpar(5)	= kappa;
  mttpar(6)	= z;
## END Code
endfunction


function [mttpar] = ChannelC_numpar;
## [mttpar] = ChannelC_numpar;
## System ChannelC, representation numpar, language m; 
## File ChannelC_numpar.m; 
## Generated by MTT on Tue Mar 10 11:29:53 AEDT 2015; 

## Horrible fudge to make mtt_m2p work
global ...
mtt_no_globals ;

   mttpar = zeros(7,1);

## BEGIN Code

## User defined code from ChannelC_numpar.txt
  c	=  1.0;
  k	=  1.0; 
  rt	=  8.3144621*300;
  f       =  96.4853399;	
  g	=  1.0; 
  kappa	=  1.0; 
  z	=  2; 

  ## Set up the mttpar vector
  mttpar(1)	= c;
  mttpar(2)	= f;
  mttpar(3)	= k;
  mttpar(4)	= rt;
  mttpar(5)	= g;
  mttpar(6)	= kappa;
  mttpar(7)	= z;
## END Code
endfunction
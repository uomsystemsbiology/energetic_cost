
function [mttpar] = GateNaC_numpar;
## [mttpar] = GateNaC_numpar;
## System GateNaC, representation numpar, language m; 
## File GateNaC_numpar.m; 
## Generated by MTT on Wed Mar 11 09:48:21 AEDT 2015; 

## Horrible fudge to make mtt_m2p work
global ...
mtt_no_globals ;

   mttpar = zeros(5,1);

## BEGIN Code

## User defined code from GateNaC_numpar.txt
  k_ch	=  1.0; 
  k_cm	=  1.0; 
  k_oh	=  1.0; 
  k_om	=  1.0; 
  rt	=  1.0; 

  ## Set up the mttpar vector
  mttpar(1)	= k_ch;
  mttpar(2)	= k_cm;
  mttpar(3)	= k_oh;
  mttpar(4)	= k_om;
  mttpar(5)	= rt;
## END Code
endfunction

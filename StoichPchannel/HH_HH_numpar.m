
function [mttpar] = HH_numpar;
## [mttpar] = HH_numpar;
## System HH, representation numpar, language m; 
## File HH_numpar.m; 
## Generated by MTT on Wed May 21 09:44:20 EST 2014; 

## Horrible fudge to make mtt_m2p work
global ...
mtt_no_globals ;

   mttpar = zeros(14,1);

## BEGIN Code

## User defined code from HH_numpar.txt
  hh =  1;
  k_k =  1/1000;
  k_n =  1/1000;
  rt =  8.3144621*300;		
  f =  96.4853399;			
  v_eq =  -65;			
  c_m	=  1;			
  g_l	=  0.3;			
  kappa_k =  0.058562;		
  kappa_n =  0.22816;		
  z	=  1;			
  if hh
    match_fac_k =  0.12143;
    match_fac_n =  3.1010;
    kappa_k =  kappa_k/match_fac_k;
    kappa_n =  kappa_n/match_fac_n;
    v_k =  -77;
    v_n =  50;
    v_l =  -54.4;
  else
    v_k =  0;
    v_n =  0;
    v_l =  -66.321;			
  endif

  ## Set up the mttpar vector
  mttpar(1)	= c_m;
  mttpar(2)	= f;
  mttpar(3)	= hh;
  mttpar(4)	= k_k;
  mttpar(5)	= k_n;
  mttpar(6)	= rt;
  mttpar(7)	= v_eq;
  mttpar(8)	= g_l;
  mttpar(9)	= kappa_k;
  mttpar(10)	= kappa_n;
  mttpar(11)	= v_k;
  mttpar(12)	= v_l;
  mttpar(13)	= v_n;
  mttpar(14)	= z;
## END Code
endfunction
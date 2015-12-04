function mttpar = ReactorTF_numpar();
% mttpar = ReactorTF_numpar();
%System ReactorTF, representation numpar, language m;
%File ReactorTF_numpar.m;
%Generated by MTT on Thu Aug 24 14:28:46 BST 2000;
%

#====== Set up the global variables ======#
global ...
     a ...
     b ...
     c ...
     c_0 ...
     c_a ...
     c_b ...
     c_p ...
     e_1 ...
     e_2 ...
     e_3 ...
     f_s ...
     h ...
     h_1 ...
     h_2 ...
     h_3 ...
     k ...
     k_1 ...
     k_2 ...
     k_3 ...
     n ...
     q ...
     q_1 ...
     q_2 ...
     q_3 ...
     q_s ...
     rho ...
     t_0 ...
     t_s ...
     v_r ...
     x1 ...
     x2 ...
     x3 ;
## Set parameters to zero
 a = 0.0;
 b = 0.0;
 c = 0.0;
 c_0 = 0.0;
 c_a = 0.0;
 c_b = 0.0;
 c_p = 0.0;
 e_1 = 0.0;
 e_2 = 0.0;
 e_3 = 0.0;
 f_s = 0.0;
 h = 0.0;
 h_1 = 0.0;
 h_2 = 0.0;
 h_3 = 0.0;
 k = 0.0;
 k_1 = 0.0;
 k_2 = 0.0;
 k_3 = 0.0;
 n = 0.0;
 q = 0.0;
 q_1 = 0.0;
 q_2 = 0.0;
 q_3 = 0.0;
 q_s = 0.0;
 rho = 0.0;
 t_0 = 0.0;
 t_s = 0.0;
 v_r = 0.0;
 x1 = 0.0;
 x2 = 0.0;
 x3 = 0.0;
 %  -*-octave-*- Put Emacs into octave-mode
 %  Numerical parameter file (ReactorTF_numpar.txt)
 %  Generated by MTT at Fri Mar  3 09:22:56 GMT 2000

 %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %  %% Version control history
 %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %  %% $Id: ReactorTF_numpar.m,v 1.1 2000/12/28 17:12:57 peterg Exp $
 %  %% $Log: ReactorTF_numpar.m,v $
 %  %% Revision 1.1  2000/12/28 17:12:57  peterg
 %  %% To RCS
 %  %%
 %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 % 
a =  0;				 %  Dummy
b =  0;				 %  Dummy
c =  0;				 %  Dummy
c_0 =  0;			 %  Dummy
c_a =  0;			 %  Dummy
c_b =  0;			 %  Dummy
c_p =  0;			 %  Dummy
e_1 =  0;			 %  Dummy
e_2 =  0;			 %  Dummy
e_3 =  0;			 %  Dummy
f_s =  0;			 %  Dummy
h =  0;				 %  Dummy
h_1 =  0;			 %  Dummy
h_2 =  0;			 %  Dummy
h_3 =  0;			 %  Dummy
k =  0;				 %  Dummy
k_1 =  0;			 %  Dummy
k_2 =  0;			 %  Dummy
k_3 =  0;			 %  Dummy
n =  0;				 %  Dummy
q =  0;				 %  Dummy
q_1 =  0;			 %  Dummy
q_2 =  0;			 %  Dummy
q_3 =  0;			 %  Dummy
q_s =  0;			 %  Dummy
rho =  0;			 %  Dummy
t_0 =  0;			 %  Dummy
t_s =  0;			 %  Dummy
v_r =  0;			 %  Dummy
x1 =  0;				 %  Dummy
x2 =  0;				 %  Dummy
x3 =  0;				 %  Dummy

 % 
rho =  900;			 %  Density
c_p =  5.0;			 %  Specific heat

 % 
k_1 =  2.5e10;			 %  Reaction rate constant
q_1 =  1e4;			 %  Exotherm constant
h_1 =  1e4;			 %  Heat of reaction

 % 
k_2 =  2.65e12;			 %  Reaction rate constant
q_2 =  1.2e4;			 %  Exotherm constant
h_2 =  1.2e4;			 %  Heat of reaction

 % 
k_3 =  6e7;			 %  Reaction rate constant
q_3 =  8e3;			 %  Exotherm constant
h_3 =  3e4;			 %  Heat of reaction

 % 
c_0 =  10;			 %  Inflow conc
t_0 =  500;			 %  Inflow temp

 % 
t_s =  530;			 %  Steady-state temp
f_s =  100;			 %  Steady-state flow








## Set up the parameter vector
  mttpar(1) 	= a;
  mttpar(2) 	= b;
  mttpar(3) 	= c;
  mttpar(4) 	= c_0;
  mttpar(5) 	= c_a;
  mttpar(6) 	= c_b;
  mttpar(7) 	= c_p;
  mttpar(8) 	= e_1;
  mttpar(9) 	= e_2;
  mttpar(10) 	= e_3;
  mttpar(11) 	= f_s;
  mttpar(12) 	= h;
  mttpar(13) 	= h_1;
  mttpar(14) 	= h_2;
  mttpar(15) 	= h_3;
  mttpar(16) 	= k;
  mttpar(17) 	= k_1;
  mttpar(18) 	= k_2;
  mttpar(19) 	= k_3;
  mttpar(20) 	= n;
  mttpar(21) 	= q;
  mttpar(22) 	= q_1;
  mttpar(23) 	= q_2;
  mttpar(24) 	= q_3;
  mttpar(25) 	= q_s;
  mttpar(26) 	= rho;
  mttpar(27) 	= t_0;
  mttpar(28) 	= t_s;
  mttpar(29) 	= v_r;
  mttpar(30) 	= x1;
  mttpar(31) 	= x2;
  mttpar(32) 	= x3;

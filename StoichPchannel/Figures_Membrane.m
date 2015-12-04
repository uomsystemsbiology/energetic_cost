Setpath;
Setplot;
graphics_toolkit("gnuplot");

## Set up state etc.
last = 20e-3;

par = Membrane_numpar;		# Parameters
sympar = Membrane_sympar;	# Symbolic parameters

Set_Membrane;


## Simulation - scaled units
[x X v v_fix Energy_var] = stoich_sim ("Membrane","",t,X_0_d, i_fix);

FV = X(1,:)./(par(sympar.C_m)/F^2); # State of electrogenic cap.
V = FV/F;			    # Voltage

## Put back actual unscaled flows
v = v_unit*v;			# mol/sec
v_fix = v_unit*v_fix;		# mol/sec
C_m = v_unit*par(sympar.C_m);

## Put back actual unscaled energy flow.
P_r = RT*v_unit*Energy_var.P_r;	# Dissipated power. W
P_e = -RT*v_unit*Energy_var.P_ce;	# External power required W
P_c = RT*v_unit*Energy_var.P_c;	# Capacitive power. W

## External power
i_ext = [1 3 5];		# Indices of external flows
P_e_k = sum(P_e(1:2,:));	# Power- k channel
P_e_n = sum(P_e(3:4,:));	# Power- Na channel
P_e_l = sum(P_e(5:6,:));	# Power- NL channel
P_e_all = sum(P_e(1:6,:));

## Internal power
i_ch = [2 5 6];
P_r_k = P_r(2,:);
P_r_n = P_r(5,:);
P_r_l = P_r(6,:);

## Chemical potentials
mu = RT*Energy_var.mu;
mu_fix = mu(i_fix,:);

## Equiv free energies (constant as concs are fixed)
G_k = mu_fix(2,1)-mu_fix(1,1);
G_n = mu_fix(4,1)-mu_fix(3,1);
G_l = mu_fix(6,1)-mu_fix(5,1);

## Affinities
A = RT*Energy_var.A;
A_k = A(2,:);
A_n = A(5,:);
A_l = A(6,:);

## Channel fluxes
v_k = v(2,:);			# mol/sec
v_n = v(5,:);
v_l = v(6,:);

## External fluxes
v_e_k = v_fix(1,:);
v_e_n = v_fix(3,:);
v_e_l = v_fix(5,:);

## Check dissipated power
check_norm(A_n.*v_n-P_r_n,"P_r_n")

## Check external power
check_norm(-v_e_n*G_n-P_e_n,"P_e_n")

## HasOttCal10 formula
v_atp = -v_n/3;			# ATP equivalent (mol)

## Cumulative fluxes
x_k = cumsum(v_k)*dt;
x_n = cumsum(v_n)*dt;
x_l = cumsum(v_l)*dt;

## HasOttCal10 formula
x_atp = -x_n/3;			# ATP equivalent (mol)
px_atp = x_atp*1e12;		# ATP equivalent (pmol)

## Channel currents 
i_k = F*v_k;			# A
i_n = F*v_n;			# A
i_l = F*v_l;			# A
i_net = i_k+i_n+i_l;

## Cumulative currents
q_k = cumsum(i_k)*dt;		# C
q_n = cumsum(i_n)*dt;
q_l = cumsum(i_l)*dt;

## CHECK
disp(sprintf("Na charge flow: %g nC/cm^2 (SenSte14 Fig 9: %i)",...
			         max(-q_n)*1e9, 1020))

## Gate currents
i_gate = [1 3 4];
I_g = F*v(i_gate,:);		# A

## Dissipated energy
E_r = cumsum(P_r,2)*dt;
E_e = cumsum(P_e,2)*dt;

## Per channel
E_e_k = cumsum(P_e_k)*dt;
E_e_n = cumsum(P_e_n)*dt;
E_e_l = cumsum(P_e_l)*dt;
E_e_all = cumsum(P_e_all)*dt;

## ATP-equivalent stuff
G_ATP = 31000;			# J/mol
P_atp = G_ATP*v_atp;		
E_atp = cumsum(P_atp)*dt;

## Check
check_norm(E_atp-x_atp*G_ATP,"E_atp")

## Reexpress for plotting
mV = V*1000;			# V in mV
mV_eq = V_eq*1000;		# V_eq in mV
mt = t*1000;			# t in msec
muv = v*1e6;			# mu mol/sec
nv = v*1e9;			# nmol/sec
muv_fix = v_fix*1e6;		# mumol/sec
nv_fix = v_fix*1e9;		# nmol/sec
muP_r = P_r*1e6;		# micro W
muP_e = P_e*1e6;		# micro W
nP_r = P_r*1e9;			# nano W
nP_e = P_e*1e9;			# nano W
nE_r = E_r*1e9;			# nJ
nE_e = E_e*1e9;			# nJ

nE_k = E_e_k*1e9;		# nJ
nE_n = E_e_n*1e9;		# nJ
nE_l = E_e_l*1e9;		# nJ
nE_all = E_e_all*1e9;		# nJ
nE_atp = E_atp*1e9;		# nJ

kmu = mu*1e-3;			# kJ/mol
kA  = A*1e-3;			# kJ/mol
nI_g = I_g*1e9;			# nA

aE_all = E_e_all/G_ATP;		# Mol of ATP
paE_all = aE_all*1e12;		# pMol of ATP

if (VaryNa==0)			# Supress in Figures_vary
  ## Plotting
  colour = 1;
  ## Action potential
  figure(10);
  one = ones(size(t));
  plot(mt,mV)
  grid;
  xlabel("t (msec)");
  ylabel("V(mV)");
  fig(Name,"V",colour);
  
  ## Channel flows
  figure(12);
  plot(mt,i_k*1e3,"-;K;",
       mt,i_n*1e3,"--;Na;");
       ##       mt,i_l*1e3,"..k;L;",
       ## mt,i_net*1e3,":;Net;");
  grid;                                 %
  xlabel("t (msec)");
  ylabel("I (mA)");
  fig(Name,"I",colour);
  
  ## Gate flows
  figure(13);
  plot(mt,nI_g(1,:),"-;n;",\
       mt,nI_g(2,:),"--;m;",
       mt,nI_g(3,:),":;h;");
  grid;
  xlabel("t (msec)");
  ylabel("I_g (nA)");
  fig(Name,"I_g",colour);
  
  I_k_0 = F*muv(2,1);
  I_n_0 = F*muv(5,1);
  I_l_0 = muv(6,1);
  
  ## Gating functions
  figure(15)
  plot(mt,X(3,:)/x_g,"-;n(t);",...
       mt,X(7,:)/x_g,"--;m(t);",...
       mt,X(8,:)/x_g,":;h(t);");
  grid;
  legend("location", "north");
  xlabel("t (msec)");
  ylabel("n(t),m(t),h(t)");
  fig(Name,"nmh",colour);
  
  
  ## External flows
  figure(20)
  plot(mt,nv_fix)
  grid;
  xlabel("t (msec)");
  ylabel("v_{fix} (nM/sec)");
  fig(Name,"v_fix",colour);
  
  ## Internal flows
  figure(21)
  plot(mt,nv(i_ch,:),mt,sum(nv(i_ch,:)))
  grid;
  xlabel("t (msec)");
  ylabel("v (nM/sec)");
  Legend = {"K" "N" "L" "Net"};
  legend(Legend);
  fig(Name,"v",colour);


## figure(22)
## plot(mt,kmu(1,:))
## grid;
## xlabel("t (msec)");
## ylabel("mu (kJ/mol)");
## fig(Name,"mu_m",colour);

  ## Affinities
  figure(23)
  plot(mt,kA(i_ch,:))
  grid;
  Legend = {"K" "N" "L" "Net"};
  legend(Legend);
  xlabel("t (msec)");
  ylabel("A (kJ/mol)");
  fig(Name,"A",colour);
  
  ## Ext and int power (total)
  figure(30)
  plot(mt,sum(muP_e),"-", mt,sum(muP_r),"--");
  grid;
  Legend = {"Ext.","Int."};
  legend(Legend);
  xlabel("t (msec)");
  ylabel("P_e,P_r (\\mu W)");
  fig(Name,"P_er",colour);
  
  ## Ext and int power (per channel)
  figure(31)
  plot(mt,P_r_n*1e6,";Na (int);",\
       mt,P_e_n*1e6,";Na (ext);",\
       mt,P_r_k*1e6,";K (int);",\
       mt,P_e_k*1e6,";K (ext);")
  grid;
  xlabel("t (msec)");
  ylabel("P_r (\\mu W)");
  fig(Name,"P_er_ch",colour);
  
  
## figure(32)
## plot(mt,sum(nE_e), mt,sum(nE_r))
## grid;
## xlabel("t (msec)");
## ylabel("E_e, E_r (nJ)");
## fig(Name,"E_er",colour);

## figure(33)
## plot(mt,nE_r(i_ch,:), mt,sum(nE_r(i_ch,:)))
## grid;
## xlabel("t (msec)");
## ylabel("E_r (nJ)");
## Legend = {"K" "N" "L" "Net"};
## legend(Legend);
## fig(Name,"E_r",colour);

  ## Total energy and proxy (nJ)
  figure(34)
  plot(mt,nE_all,"-;Actual;",  mt, nE_atp,"--;Proxy;")
  grid;
  legend("location","southeast")
  xlabel("t (msec)");
  ylabel("Energy (nJ)");
  fig(Name,"ATP",colour);
  
  ## Total energy and proxy (mol)
  figure(35)
  plot(mt,paE_all,";Actual;",  mt, px_atp,";Proxy;")
  grid;
  legend("location","southeast")
  xlabel("t (msec)");
  ylabel("Energy (pmol ATP)");
  fig(Name,"ATP_mol",colour);

## Plot electical energy in capacitor
 figure(36)
 V_e = V(length(mt));
 plot(mt, (0.5*C_m*(V-V_eq).^2)*1e9 );
 grid;
 xlabel("t (msec)");
 ylabel("E_c (nJ)");
 fig(Name,"E_c",colour);
endif

## Energy calculations
maxE_e = max(nE_all);
maxE_e_n = max(nE_n);
maxE_e_n = max(nE_n);
maxE_e_k = max(nE_k);
maxE_atp = max(nE_atp);
disc = 100*(maxE_atp-maxE_e)/maxE_e;
disp(sprintf("E_e=%g nJ, E_atp=%g nJ, discrepancy=%gpc",maxE_e,maxE_atp,disc));

pX_k = max(abs(x_k))*1e12
pX_n = max(abs(x_n))*1e12
pX_l = max(abs(x_l))*1e12
pX_atp = max(abs(x_atp))*1e12

disp(sprintf("x_k=%g, x_n=%g, x_l=%g",pX_k,pX_n,pX_l));
nE_comp_k = abs(G_k)*pX_k/1000;
check_norm(nE_comp_k-maxE_e_k,"E_k");

nE_comp_n = abs(G_n)*pX_n/1000;
check_norm(nE_comp_n-maxE_e_n,"E_n");

nE_comp = nE_comp_k+nE_comp_n;
check_norm((nE_comp-maxE_e)./maxE_e,"E",1e-4);

disp(sprintf("Computed: E_e=%g nJ, E_k=%g nJ, E_n=%g nJ", ...
	     nE_comp, nE_comp_k, nE_comp_n));

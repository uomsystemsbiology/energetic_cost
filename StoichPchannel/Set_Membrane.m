## [rtf R F_0] = RTF;
## V_N = RTF*1000;		# Nernst voltage in mv
## Get some constants
if !exist("VaryNa")
   VaryNa=0;
endif

const = ThermoConstants;
F = const.F;
V_N = const.V_N;
Avogadro = const.A;


## Initial depolarisation
depolarisation = 20e-3;		# V


t = [0:last/1000:last];		# Time (sec)
dt = mean(diff(t));	# Time step (sec);
i_fix = [4;5;10;11;12;13];	# Simulation
i_fix_0 = [1;4;5;10;11];	# Equilibrium
i_fix_1 = [1;4;5;10;11;12;13];	# Steady-state with depolarisation.


## Volumes
V_i = par(sympar.V_i);
V_e = par(sympar.V_e);

## Thermo. gains
K_k = par(sympar.K_k);
K_n = par(sympar.K_n);
K_l = par(sympar.K_l);

V_eq_0 = par(sympar.V_eq);
C_m = par(sympar.C_m);
RT = par(sympar.RT);
F = par(sympar.F);
kappa_k = par(sympar.kappa_k);
kappa_n = par(sympar.kappa_n);
kappa_l = par(sympar.kappa_n);
v_k = par(sympar.v_k);
v_n = par(sympar.v_n);
v_l = par(sympar.v_l);

x_g = par(sympar.x_g);		# Total gate state

v_unit = par(sympar.v_unit);	# Scaling factor

Name = "Memb";

## ## Stoichiometry
## stoich = dm2stoich("Membrane","");

## Initial conditions
x_i_k = V_i*397/K_k;		# mmol
x_e_k = V_e*20/K_k;		# mmol

x_i_n = V_i*50/K_n;		# mmol
if VaryNa
   x_i_n = vary*x_i_n;
endif
x_e_n = V_e*437/K_n;		# mmol

x_i_l = V_i*100e-3/K_l;
x_e_l = V_e*12.193e-3/K_l;		# From Figures_GHK

X_00 = [(V_eq_0)*C_m/F;0;x_g;x_i_k; x_e_k;0;x_g;0;x_g;x_i_n; x_e_n;x_i_l;x_e_l];

## Find initial steady state
[x_0 X_0 v_0] = stoich_sim ("Membrane","",0,X_00,i_fix_0);
V_eq = F*(X_0(1)/C_m);
n_0 = X_0(3);
m_0 = X_0(7);
h_0 = X_0(8);
x_i_l_0 = X_0(12);
x_e_l_0 = X_0(13);
V_l_0 = V_N*log(x_e_l_0/x_i_l_0);

## Put in depolarisation
X_0_d = X_0;
X_0_d(1) = X_0_d(1) + depolarisation*(C_m/F);

## [x_1 X_1 v_1] = stoich_sim ("Membrane","",0,X_0_d,i_fix_1);
## V_1 = F*(X_1(1)/C_m)
## n_1 = X_1(3)
## m_1 = X_1(7)
## h_1 = X_1(8)

## The Free-energy constants
K_c = Membrane_K(par);

## Compute some numbers
K_m = (F^2)/(C_m*RT);


Setpath;
Setplot;
graphics_toolkit("gnuplot");

## Set up state etc.
last = 20e-3;

par = Membrane_numpar;		# Parameters
sympar = Membrane_sympar;	# Symbolic parameters

Set_Membrane;


## Simulation

[x X v v_fix Energy_var] = stoich_sim ("Membrane","",t,X_0_d, i_fix);

FV = X(1,:)./(par(sympar.C_m)/F^2); # State of electrogenic cap.
V = FV/F;			    # Voltage

## Put back actual unscaled flowsp
v = v_unit*v;
v_fix = v_unit*v_fix;

## Put back unscaled energy flow.
P_r = v_unit*Energy_var.P_r;	# Dissipated power.
P_e = -v_unit*Energy_var.P_ce;	# External power required

## Channel fluxes
i_ch = [2 5 6];
v_k = v(2,:);
v_n = v(5,:);
v_l = v(6,:);


## Gate currents
i_gate = [1 3 4];
I_g = F*v(i_gate,:);

## Dissipated energy
dt = mean(diff(t));	# Sec.
E_r = cumsum(P_r,2)*dt;
E_e = cumsum(P_e,2)*dt;

## ATP stuff
G_ATP = 31000;			# J/mol
E_ATP = G_ATP/Avogadro;		# J per molecule
v_ATP = -v_n/3;			# Assume 3 Na per ATP
x_ATP = cumsum(v_ATP)*dt;	# moles
N_ATP = x_ATP*Avogadro;
GN_ATP = N_ATP/1e9		# G ATPs

## Reexpress for plotting
mV = V*1000;			# V in mV
mV_eq = V_eq*1000;		# V_eq in mV
mt = t *1000;			# t in msec
muv = v*1e6;			# mu mol/sec
nv = v*1e9;			# nmol/sec
muv_fix = v_fix*1e6;		# mumol/sec
nv_fix = v_fix*1e9;		# nmol/sec
muP_r = P_r*1e6;		# micro W
muP_e = P_e*1e6;		# micro W
nP_r = P_r*1e9;			# nano W
nP_e = P_e*1e9;			# nano W
pE_r = E_r*1e12;		# pJ
pE_e = E_e*1e12;		# pJ
aE_r = E_r/E_ATP;		# Energy per mol of ATP
aE_e = E_e/E_ATP;		# Energy per mol of ATP
GaE_r = aE_r/1e9;		# Energy per mol of ATP Giga
GaE_e = aE_e/1e9;		# Energy per mol of ATP Giga
kmu = Energy_var.mu*1e-3;	# kJ/mol
kA  = Energy_var.A*1e-3;	# kJ/mol
muI_g = I_g*1e6;		# muA

## Specific channels
nP_K = sum(nP_e([1:2],:));	# nW
nP_N = sum(nP_e([3:4],:));	# nW
nP_L = sum(nP_e([5:6],:));	# nW

pE_K = sum(pE_e([1:2],:));	# pJ
pE_N = sum(pE_e([3:4],:));	# pJ

pA_K = pE_K/G_ATP;		# nmol ATP
pA_N = pE_N/G_ATP;		# nmol ATP




## Channel currents (mu A)
I_k = F*muv(2,:);
I_n = F*muv(5,:);
I_l = F*muv(6,:);
I = I_k + I_n + I_l;

figure(10);
one = ones(size(t));
plot(mt,mV,mt,one*mV_eq)
grid;
xlabel("t (msec)");
ylabel("V(mV)");
##axis([0 2 -100 150])
fig(Name,"V",2);

## Work out cumulative external flows
Total_mol = sum(v_fix([1:2:6],:)')*dt # Moles
Total = Total_mol*Avogadro

## Cumulative disipated power
Energy_r = sum(sum(P_r(i_ch,:)))*dt # Joules
nEnergy_r = Energy_r*1e9		  # nJ



figure(11);
plot(mt,I);
grid;
xlabel("t (msec)");
ylabel("I_{net} (\\mu A)");
fig(Name,"I_all",2);

figure(12);
plot(mt,I_k,";K;",
     mt,I_n,";Na;",
     mt,I_l,";L;",
     mt,I,";Net;");
grid;
xlabel("t (msec)");
ylabel("I (\\mu A)");
fig(Name,"I",2);

figure(13);
plot(mt,muI_g(1,:),";n;",\
     mt,muI_g(2,:),";m;",
     mt,muI_g(3,:),";h;");
grid;
xlabel("t (msec)");
ylabel("I_g (\\mu A)");
fig(Name,"I_g",2);

I_k_0 = F*muv(2,1)
I_n_0 = F*muv(5,1)
I_l_0 = muv(6,1)

figure(15)
plot(mt,X(3,:)/x_g,";n(t);",...
     mt,X(7,:)/x_g,";m(t);",...
     mt,X(8,:)/x_g,";h(t);");
grid;
xlabel("t (msec)");
ylabel("n(t),m(t),h(t)");
fig(Name,"nmh",2);


## ## Flows are in M/sec
## ## Mu is in J/mol
figure(20)
plot(mt,nv_fix)
grid;
xlabel("t (msec)");
ylabel("v_{fix} (nM/sec)");
fig(Name,"v_fix",2);

figure(21)
plot(mt,nv(i_ch,:),mt,sum(nv(i_ch,:)))
grid;
xlabel("t (msec)");
ylabel("v (nM/sec)");
Legend = {"K" "N" "L" "Net"};
legend(Legend);
fig(Name,"v",2);

figure(22)
plot(mt,kmu(1,:))
grid;
xlabel("t (msec)");
ylabel("mu (kJ/mol)");
fig(Name,"mu_m",2);

figure(23)
plot(mt,kA(i_ch,:))
grid;
Legend = {"K" "N" "L" "Net"};
legend(Legend);
xlabel("t (msec)");
ylabel("A (kJ/mol)");
fig(Name,"A",2);

figure(30)
plot(mt,sum(muP_e),mt,sum(muP_r));
grid;
Legend = {"Ext.","Chan."};
legend(Legend);
xlabel("t (msec)");
ylabel("P_e,P_r (\\mu W)");
fig(Name,"P_e",2);

figure(31)
plot(mt,muP_r(i_ch,:), mt,sum(muP_r(i_ch,:)))
grid;
Legend = {"K" "N" "L" "Net"};
legend(Legend);
xlabel("t (msec)");
ylabel("P_r (\\mu W)");
fig(Name,"P_r",2);

figure(32)
plot(mt,sum(pE_e), mt,sum(pE_r))
grid;
xlabel("t (msec)");
ylabel("E_e, E_r (pJ)");
fig(Name,"E_er",2);

figure(33)
plot(mt,pE_r(i_ch,:), mt,sum(pE_r(i_ch,:)))
grid;
xlabel("t (msec)");
ylabel("E_r (pJ)");
Legend = {"K" "N" "L" "Net"};
legend(Legend);
fig(Name,"E_r",2);

figure(34)
plot(mt,sum(GaE_e),";Actual;",  mt, GN_ATP,";Proxy;")
grid;
axis([0,max(mt),0,5e9]);
xlabel("t (msec)");
ylabel("ATP molecules");
fig(Name,"ATP",2);


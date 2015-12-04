Setpath; Setplot


const = ThermoConstants;
F = const.F;

## Parameters etc.
par = ChannelNaC_numpar;
sym = ChannelNaC_sympar;
K_n = par(sym.K_n);
v_n = 115e-3;
x_g = par(sym.x_g)			# Total gate state

V = [-100:-41 -39:100]*1e-3;	# Actual voltage
V_n = V-v_n;			# Relative voltage

## Constant states
x_i_n = 50/K_n;			# mM
x_e_n = 437/K_n;

X_0 = [F*V(1)
       x_g
       0
       x_g
       0
       x_i_n
       x_e_n];	
 
i_fix = [1;6;7];

m = [];
h = [];
i_n = [];
for i = 1:length(V)
  X_0(1) = F*V(i);
  [x X v] = stoich_sim ("ChannelNaC","",0,X_0,i_fix);
  x_m = X(3);
  x_h = X(4);
  m = [m x_m/x_g];
  h = [h x_h/x_g];
  i_n = [i_n F*v(3)];
end


## HH equivalent
g_hh = 36e-3;
i_hh = g_hh * V_n;
i_hh_gated = (m.^3).*h.*i_hh;

## Scale
mV = 1000*V;			# mV
mi_n = 1000*i_n;		# mA
mi_hh_gated = 1000*i_hh_gated;	# mA

mu = F*V;
v_n  = i_n/F;
v_hh  = i_hh_gated/F;

kmu = mu*1e-3;			# kJ/mol
nv_n = v_n*1e9;			# nmol/sec
nv_hh = v_hh*1e9;		# nmol/sec

figure(20);
plot(kmu,m, ";m;", kmu,h, "--;h;" );
xlabel("\\mu (kJ/mol)");
ylabel("m,h");
grid
fig("Mn","vx_ss",2)

figure(21);
plot(kmu,nv_n,";GHK;", kmu,nv_hh,"--;HH;")
xlabel("\\mu (kJ/mol)");
ylabel("v_n, v_{hh} (nmol/sec)");
grid
fig("Mn","vi_ss",2)

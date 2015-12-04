global U

const = ThermoConstants;
F = const.F;

## Parameters etc.
par = ChannelKC_numpar;
sym = ChannelKC_sympar;
K_k = par(sym.K_k)
v_k = -77e-3;
x_g = par(sym.x_g)			# Total gate state

##V_k = [-50:1:150];		# Relative voltage
##V = V_k + v_k;			# Actual voltage
V = [-100:-56 -54:100]*1e-3;		# Miss -55 ????
V_k = V-v_k;

## Constant states
x_i_k = 397/K_k;			# mM
x_e_k = 20/K_k;	

X_0 = [F*V(1)
       x_g
       0
       x_i_k
       x_e_k];

i_fix = [1;4;5];
 
n_ss = [];
i_k = [];
for i = 1:length(V)
  X_0(1) = F*V(i);
  [x X v] = stoich_sim ("ChannelKC","",0,X_0,i_fix);
  i_k_i = F*(v(1) +v(2));
  n_ss = [n_ss, X(3)/x_g];
  i_k = [i_k, i_k_i];
end


## HH equivalent
g_hh = 36e-3;
i_hh = g_hh * V_k;
i_hh_gated = (n_ss.^4).*i_hh;

mV = 1000*V;			# mV
mi_k = 1000*i_k;		# mA
mi_hh_gated = 1000*i_hh_gated;	# mA

mu = F*V;
v_k  = i_k/F;
v_hh  = i_hh_gated/F;

kmu = mu*1e-3;			# kJ/mol
nv_k = v_k*1e9;			# nmol/sec
nv_hh = v_hh*1e9;		# nmol/sec


figure(20);
plot(kmu,n_ss);
xlabel("\\mu (kJ/mol)");
ylabel("n");
grid
fig("Mk","vx_ss",2)

figure(21);
plot(kmu,nv_k,";GHK;" , kmu,nv_hh,";HH;")
xlabel("\\mu (kJ/mol)");
ylabel("v_k, v_{hh} (nmol/sec)");
grid
fig("Mk","vi_ss",2)

## figure(22);
## plot(V_k,i_k./V_k)
## xlabel("v_k")
## ylabel("g");
## grid
## fig("Mk","g",2)
## ## 

const = ThermoConstants;
F = const.F;

## Parameters etc.
par = ChannelKC_numpar;
sym = ChannelKC_sympar;
K_k = par(sym.K_k)
v_k = -77;

x_g = par(sym.x_g)			# Total gate state

## Initial cond.
## Constant states
x_i_k = 397/K_k;			# mM
x_e_k = 20/K_k;	

X_0 = [0
       x_g
       0
       x_i_k
       x_e_k];


## simulate states
i_fix = [1;4;5];
t = [0:0.1:50]'*1e-3;		# Time
[x X v] = stoich_sim ("ChannelKC","",t,X_0,i_fix,[],"","ChannelKC_Xfix");
i_k = F*v(2,:);
one = ones(size(t));
V = 0.1*sin(2*pi*1e3*t/10);	# V
x_open = X(3,:);

mV = 1000*V;			# mV
mi_k = 1000*i_k;		# mV

mu = F*V;
v_k  = i_k/F;
kmu = mu*1e-3;			# kJ/mol
nv_k = v_k*1e9;			# nmol/sec

## Plot
R = find(t>max(t)/2);
figure(10);
plot(kmu(R),nv_k(R));
grid;
xlabel("\\mu (kJ/mol)");
ylabel("v_k (nmol/sec)");
fig("Mk","vi_sin",2)

figure(11);
plot(kmu(R),x_open(R)/x_g);
grid;
xlabel("\\mu (kJ/mol)");
ylabel("n");
fig("Mk","vx_sin",2)

## figure(12);
## plot(t,i_k);
## grid;
## xlabel("t");
## ylabel("i_k");
## fig("Mk","i_sin",2)



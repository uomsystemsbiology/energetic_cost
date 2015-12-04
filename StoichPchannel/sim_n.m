const = ThermoConstants;
F = const.F;

## Parameters etc.
par = ChannelNaC_numpar;
sym = ChannelNaC_sympar;
K_n = par(sym.K_n)
x_g = par(sym.x_g)			# Total gate state

## Initial states
x_i_n = 50/K_n;			# mM
x_e_n = 437/K_n;

X_0 = [F*0
       x_g
       0
       x_g
       0
       x_i_n
       x_e_n];	
 
i_fix = [1;6;7];

## Simulate
t = [0:0.1:50]*1e-3;
[x X v] = stoich_sim ("ChannelNaC","",t,X_0,i_fix,[],"","ChannelKC_Xfix");

m = X(3,:)/x_g;
h = X(4,:)/x_g;
i_n = F*v(3,:);
V = 0.1*sin(2*pi*1e3*t/10);		# V

## Plot
R = find(t>max(t)/2);

## v_n = 115;
## i_n = Y(1,:);
## V = U;
mV = 1000*V;			# mV
mi_n = 1000*i_n;		# mA

mu = F*V;
v_n  = i_n/F;
kmu = mu*1e-3;			# kJ/mol
nv_n = v_n*1e9;			# nmol/sec

figure(10);
plot(kmu(R),nv_n(R));
grid;
xlabel("\\mu (kJ/mol)");
ylabel("v_n (nmol/sec)");
fig("Mn","vi_sin",2)

figure(11);
plot(kmu(R),m(R),";m;", kmu(R),h(R),"--;h;");
grid;
xlabel("\\mu (kJ/mol)");
ylabel("m,h");
fig("Mn","vx_sin",2)


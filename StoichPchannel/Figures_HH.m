Setpath;
Setplot;

[rtf R F_0] = RTF;
V_N = RTF*1000;		# Nernst voltage in mv

last = 50;
depolarisation = 7;		# mV
depolarisation = 15;		# mV
##  depolarisation = 50;		# mV
depolarisation = 60;		# mV
##depolarisation = 10;		# mV
## depolarisation = 80;		# mV

t = [0:0.01:last];
i_fix = [4;5;10;11];

par = HH_numpar;		# Parameters
sympar = HH_sympar;		# Symbolic parameters

##HH = par(sympar.HH)

K_k = par(sympar.K_k)
K_n = par(sympar.K_n)

V_eq_0 = par(sympar.V_eq)
C_m = par(sympar.C_m);
RT = par(sympar.RT);
F = par(sympar.F);
kappa_k = par(sympar.kappa_k);
kappa_n = par(sympar.kappa_n);
g_l = par(sympar.g_l);
v_k = par(sympar.v_k);
v_n = par(sympar.v_n);
v_l = par(sympar.v_l)

HH=0;
if HH
  name = "GHK_HH";
else
  name = "GHK";
endif

Name = sprintf("%s",name);

## Initial conditions
x_i_k = 397/K_k;		# mM
x_e_k = 20/K_k;			# mM

x_i_n = 50/K_n;			# mM
x_e_n = 437/K_n;		# mM

## [dx gate_K_x_0] = mk_cr (0,V_eq);
## [dx gate_Na_1_x_0] = mn_cr1 (0,V_eq)
## [dx gate_Na_2_x_0] = mn_cr2 (0,V_eq)

## GKH = ghk_fun(V_eq)

## i_k_0 = GKH*F*kappa_k*gate_K_x_0^4*K_k*(exp(V_eq/V_N)*x_i_k - x_e_k)
## i_n_0 = GKH*F*kappa_n*gate_Na_1_x_0^3*gate_Na_2_x_0*K_n*(exp(V_eq/V_N)*x_i_n - x_e_n)
## i_l_0 = -(i_k_0 + i_n_0)		# Equilibrium leakage

## V_l_0 = V_eq - i_l_0/g_l

## CHECK_v_l = norm(V_l_0-v_l)

## X_00_alt = [(V_eq_0)*C_m/F
##        (1-gate_K_x_0); gate_K_x_0; 
##        x_i_k; x_e_k; 
##        (1-gate_Na_1_x_0); gate_Na_1_x_0; 
##        gate_Na_2_x_0; (1-gate_Na_2_x_0);
##        x_i_n; x_e_n];

X_000 = [(V_eq_0)*C_m/F;0;1;x_i_k; x_e_k;0;1;0;1;x_i_n; x_e_n];

## Find initial steady state
i_fix_0 = [1;i_fix];
[x_00 X_00 v_00] = stoich_sim ("HH","",0,X_000,i_fix_0);
V_eq = F*(X_00(1)/C_m)

##CHECK_X_00 = norm(X_00 - X_00_alt)

## Put in depolarisation
X_0 = X_00;
X_0(1) = X_0(1) + depolarisation*(C_m/F);

## X = lsode("HH_fun",X_0,t)';

## Y = [];
## for i = 1:length(t)
##   xx = zeros(4,1);			# Not needed
##   yy = 0;			# Not needed
##   y = HH_odeo(X(:,i),yy,t(i),par);
##   Y = [Y y];
## end

## V = Y(1,:);
[x X v] = stoich_sim ("HH","",t,X_0, i_fix);

FV = X(1,:)./(par(sympar.C_m)/F^2);
V = FV/F;
figure(10);
plot(t,V-V_eq)
grid;
xlabel("t (msec)");
ylabel("V-V_{eq}");
##axis([0 2 -100 150])
fig(Name,"V",2);

# figure(11);
# plot(t,exp(Y([3,5],:)/RT))
# grid;
# xlabel("t (msec)");
# ylabel("Gate");

## Channel currents 
I_k = F*v(2,:);
I_n = F*v(5,:);
I_l = F*v(6,:);
I = I_k + I_n + I_l;

## Gate currents
I_g = F*v([1 3 4],:);


figure(11);
plot(t,I);
grid;
xlabel("t (msec)");
ylabel("I_{net}");
fig(Name,"I_all",2);

figure(12);
plot(t,I_k,";K;",\
     t,I_n,";Na;",
     t,I_l,";L;");
grid;
xlabel("t (msec)");
ylabel("I");
fig(Name,"I",2);

figure(13);
plot(t,I_g(1,:),";m;",\
     t,I_g(2,:),";n;",
     t,I_g(3,:),";h;");
grid;
xlabel("t (msec)");
ylabel("I_g");
fig(Name,"I",2);

I_k_0 = F*v(2,1)
I_n_0 = F*v(5,1)
I_l_0 = v(6,1)

## if HH
##  V_k = v_k;
##  V_n = v_n;
## else
##  V_k = -77.251;
##  V_n = 56.045;
## endif

## g_k = I_k./(V-V_k);
## g_n = I_n./(V-V_n);
## one = ones(size(t));

## figure(14)
## plot(t,g_k,";K;",\
##      t,g_n,";Na;",\
##      t,g_l*one,";L;");
## grid;
## xlabel("t (msec)");
## ylabel("g");
## fig(Name,"g",2);

figure(15)
plot(t,X(3,:),";n(t);",\
     t,X(7,:),";m(t);",\
     t,X(8,:),";h(t);");
grid;
xlabel("t (msec)");
ylabel("n(t),m(t),h(t)");
fig(Name,"nmh",2);

## figure(16)
## plot(t,X(5,:)*K_k,";c_{ek};",\
##      t,X(11,:)*K_n,";c_{in};");
## grid;
## xlabel("t (msec)");
## ylabel("x");
## fig(Name,"x",2);


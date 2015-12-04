## Plot the GHK function
Setplot;
Setpath;

if !exist("Ion")
  Ion = "K";
endif
Ion

V_eq = -65e-3;			# mV from K&S

## Get some constants
const = ThermoConstants;
F = const.F;
V_N = const.V_N;

if strcmp(Ion,"Na")
  ## Values from KeeSne09 Table 2.1
  c_i = 50;			# mM
  c_e = 437;			# mM
  g = 120e-3;			# S/cm^2
  ## Guess
  k_mm = 0.1;
 elseif strcmp(Ion,"K")
  ## Values from KeeSne09 Table 2.1
  c_i = 397;
  c_e = 20;
  g = 36e-3;			# S/cm^2
  ## Guess
  k_mm = 50;
 elseif strcmp(Ion,"L")		# Leakage
   c_i = 100;
   v_0 = -54.4e-3;		# Give value
   c_e = c_i*exp(v_0/V_N)	# Equiv. concentration
   g = 0.3e-3;			# S/cm^2
   ## Guess
   k_mm = 50;
 else
   error(sprintf("Ion %s not recognised", Ion));
 endif

V_0 = V_N*log(c_e/c_i)
V_match = -V_0;		# Voltage at which HH & MM matches GHK
FV_match = F*V_match;

V = [-99:1:99]*1e-3;		# V
V = V + 1e-6;			# Avoid zero!
VV = V/V_N;			# V_bar = V/V_N

FV = F*V;

## Compute HH flow
VT = (V-V_0)/V_N;		# V tilde
VT_match = (V_match-V_0)/V_N;	# V tilde at match voltage

v_HH = g*(V_N/F)*VT;		# HH
v_HH_match =  g*(V_N/F)*VT_match; # Matched

## Work out equivalent kappa. (Match at V_match)
## BG
v_BG_match = (exp(VT_match) - 1); 
kappa_bg = (v_HH_match/v_BG_match)*(1/c_e);

## GHK
kappa_ghk = kappa_bg/ghk_fun(F*V_match,0)

## MM
kappa_mm = kappa_bg/mm_fun(F*V_match,V_0,k_mm)

## HH Using BG
kappa_hh = kappa_bg/ghk_fun(F*V_match,V_0)

## Compute flows.
v_0 =   (exp(VT) - 1);
v_BG =  kappa_bg*c_e*v_0;
v_GHK = kappa_ghk*c_e*ghk_fun(F*V,0).*v_0;
v_MM  = kappa_mm*c_e*mm_fun(F*V,V_0,k_mm).*v_0;

## Check
v_BG_HH = kappa_hh*c_e*ghk_fun(F*V,V_0).*v_0;



## Plot the flows.
colour = 1;
figure(1)
plot(V*1e3,v_GHK*1e9, "-;GHK;", ...
     V*1e3,v_BG_HH*1e9, "--;HH;");
grid;
Xlabel = "V (mV)";
title(Ion)
xlabel(Xlabel)
ylabel("v (nmol/sec)")
legend("location","northwest")
fig(Ion,"v",colour)

## Plot the currents
figure(2)
plot(V*1e3,F*v_GHK*1e3, "-;GHK;", ...
     V*1e3,F*v_BG_HH*1e3, "--;HH;");
grid;
Xlabel = "V (mV)";
title(Ion)
xlabel(Xlabel)
ylabel("i mA")
legend("location","northwest")
fig(Ion,"i",colour)

## Plot the Ratios.
figure(3)
plot(1000*V,v_GHK./v_BG, "-;GHK;", ...
     1000*V,v_BG_HH./v_BG, "--;HH;");
grid;
title(Ion)
xlabel(Xlabel)
ylabel("G")
legend("location","northwest")
fig(Ion,"G",colour)


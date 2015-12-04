Setpath;
Setplot;
graphics_toolkit("gnuplot");

## Set up state etc.
last = 20e-3;
Set_Membrane;
Name = "Memb_no"

## Don't fix the concentrations
[x X v v_fix Energy_var] = stoich_sim ("Membrane","",t,X_0_d);

FV = X(1,:)./(par(sympar.C_m)/F^2); # State of electrogenic cap.
V = FV/F;			    # Voltage

## Put back actual unscaled flows
v = v_unit*v;
v_fix = v_unit*v_fix;

## Put back unscaled energy flow.
P_r = v_unit*Energy_var.P_r;	# Dissipated power.
P_e = v_unit*Energy_var.P_ce;	# External power.

## Reexpress for plotting
mV = V*1000;			# V in mV
mV_eq = V_eq*1000;		# V_eq in mV
mt = t *1000;			# t in msec
muv = v*1e6;			# mu mol/sec
nv = v*1e9;			# nmol/sec

## Concentrations
i_conc = i_fix;
i_i = i_conc(1:2:6);		# Internal
i_e = i_conc(2:2:6);		# External

conc_i = X(i_i,:)/V_i;
conc_e = X(i_e,:)/V_e;


figure(10);
one = ones(size(t));
plot(mt,mV,mt,one*mV_eq)
grid;
xlabel("t (msec)");
ylabel("V(mV)");
fig(Name,"V",2);

figure(11)
plot(mt,conc_e(1:2,:))



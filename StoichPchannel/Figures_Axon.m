Setpath;
Setplot;
graphics_toolkit("gnuplot");

## Set up state etc.
last = 25e-3;;
par = Axon_numpar;		# Parameters
sympar = Axon_sympar;	# Symbolic parameters

Set_Membrane;
Name = "Axon";

## Stoichiometry
stoich_memb = dm2stoich("Membrane","");
stoich_axon = dm2stoich("Axon","");

## Sizes
N_X = stoich_axon.n_X;
n_X = stoich_memb.n_X;
N_seg = N_X/n_X

n_v = stoich_memb.n_V + 1;	# Extra R
N_v = N_seg*n_v - 1;		

XX_00 = [];
ii_fix_0 = [];
for i=1:N_seg
    XX_00 = [XX_00; X_0];
    ii_fix_0 = [ii_fix_0; (i_fix_0 + (i-1)*n_X)];
endfor


## Refine  XX_0
##[xx_0 XX_0 vv_0] = stoich_sim ("Axon","",0,XX_00,ii_fix_0);
XX_0 = XX_00;

## Put in depolarisation
XX_0_d = XX_0;
dep = depolarisation*(C_m/F);
XX_0_d(1) = XX_0_d(1) + dep;

## And add a small amount
 ## I_1 = 1:n_X:N_X
 ## XX_0_d(I_1) = XX_0_d(I_1) + 0.2*dep*ones(N_seg,1);

## Simulation
[x X v] = stoich_sim ("Axon","",t,XX_0_d, []);

FV = X(1:n_X:N_X,:)./(par(sympar.C_m)/F^2);
V = FV/F;

## Figure scaling
mt = t*1e3;			# msec
mV = V*1e3;			# mV
mV_eq = V_eq*1e3;		# mV

figure(10);
one = ones(size(t));
plot(mt,mV, mt,mV_eq*one,".")
grid;
xlabel("t (msec)");
ylabel("V (mV)");
Legend = num2str([1:N_seg]')
legend(Legend)
##axis([0 2 -100 150])
fig(Name,"V",2);

figure(11);
i_seg = 1:N_seg;
mesh(mt,i_seg,mV)
xlabel("t (msec)");
ylabel("i_{seg}");
zlabel("V (mV)");

fig(Name,"V3D",2);


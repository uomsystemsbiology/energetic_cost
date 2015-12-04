## Compare simple physical 2-state model (KeeSne 3.5.1) with HH model.
Setplot; Setpath
graphics_toolkit('gnuplot')
addpath("~/WORK/Research/SystemsBiology/Notes/2015/IonChannels/Examples/StoichPchannel");
i_fig = 10;

## Get some constants
const = ThermoConstants;
F = const.F;
V_N = const.V_N;

## KeeSne 5.1 subtract V_eq in their equations
V_eq = -65; 			# mV (KeeSne 5.1  p 206)
V_0 = 0;
v = [-10:0.1 0.1:0.1:100];	# KeeSne Voltage (mV)
V = (v + V_eq)*1e-3;;			# actual voltage (V)
one = ones(size(v));

##K = 100;			# C gain 
K = 1;				# C gain ?????
Gate_types = {"n" "m" "h"};

for i = 1:length(Gate_types);
  gate_type = Gate_types{i};
  disp(" ")
  disp("###############")
  name = sprintf("Pgate_%s", gate_type)
  disp("###############")

  ## Hodgkin-Huxley equations + Physical parameters.
  if strcmp(gate_type,"n")
    V_g = V_N;

    HH.alpha = 0.01*(-v+10)./(exp( (-v+10)/10 ) - 1);
    HH.beta  = 0.125*exp(-v/80); 

  elseif strcmp(gate_type,"m")
    V_g = V_N/3;

    HH.alpha = 0.1*(-v+25)./(exp( (-v+25)/10 ) - 1);
    HH.beta  = 4*exp(-v/18);

  elseif strcmp(gate_type,"h")
    V_g = V_N/4;

    HH.beta = 0.07*exp(-v/20);
    HH.alpha  = 1./(exp((30-v)/10) + 1);
  else
      error(sprintf("Gate type %s not known",gate_type))
  endif

  ## HH stuff in msec - change to sec.
  HH.alpha = HH.alpha*1e3;
  HH.beta = HH.beta*1e3;

  ## Gain and time-constant
  [HH.g,HH.tau] = GateProps(HH.alpha,HH.beta);
  
  ## Physical (unit parameters)
  ph.alpha = exp(V/V_g);
  ph.beta = one;
  
  ## Rescale
  i_0 = find(v==V_0);
  factor = (HH.beta(i_0)/HH.alpha(i_0))/(ph.beta(i_0)/ph.alpha(i_0));
  ph.beta =  ph.beta*factor;

  ## Make k for gate output unity.
  if strcmp(gate_type,"h")
    k_c = K;
    k_o = K*factor;
  else
    k_c = K/factor;
    k_o = K;
  endif


  k_c = k_c
  k_o = k_o

  K_eq = k_c/k_o;
  
  ## Recompute
  ph.alpha = k_c*exp(V/V_g);
  ph.beta  = k_o*one;

  ## Gain and time-constant
  [ph.g,ph.tau] = GateProps(ph.alpha,ph.beta);
  
  ## Work out kappa function
  kappa = ph.tau./HH.tau;	# Ideal

  ## CHECK
  kappa_mma = mma_cr(i,F*V);
  CHECK_mma_cr = max(abs(kappa-kappa_mma))/ max(abs(kappa_mma))

  ## Multiply and recompute
  [gph.g,gph.tau] = GateProps(kappa.*ph.alpha, kappa.*ph.beta);
  
  ##Plots
  colour = 1;
  figure(i_fig++); 
  if strcmp(gate_type,"h")
    plot(V*1e3,1-gph.g,"-;Phy;", V*1e3,1-HH.g,"--;HH;");
  else
    plot(V*1e3,gph.g,"-;Phy;", V*1e3,HH.g,"--;HH;");
  endif
  y_lim = ylim;
  axis([min(V*1e3) max(V*1e3) 0 y_lim(2)])
  xlabel("V (mV)")
  ylabel("g_{ss}")
  title(gate_type);
  legend("location","southeast")
  grid;
  fig(name, "g",colour);
  
  figure(i_fig++); 
  plot(V*1e3,gph.tau*1e3,"-;Phy;", V*1e3,HH.tau*1e3,";--HH;");
  y_lim = ylim;
  axis([min(V*1e3) max(V*1e3) 0 y_lim(2)])
  xlabel("V (mV)")
  ylabel("\\tau (msec)")
  title(gate_type);
  grid;
  fig(name, "tau",colour);
  
  ## figure(12); 
  ## plot(V,kappa);
  ## xlabel("V")
  ## ylabel("\\kappa")
  ## grid
  ## y_lim = ylim;
  ## axis([min(V) max(V) 0 y_lim(2)])
  ## fig(name, "kappa",colour);
  
endfor


clear;
Setpath;
Setplot;
graphics_toolkit("gnuplot")

## Channel design
Ions = {"K" "Na" "L"};
for i=1:length(Ions)
    Ion = Ions{i}
    Figures_GHK;
endfor

## Gate design
CompareModel;

## Memristor stuff
Figures_n;
Figures_k;

## Membrane simulation - vary Na
clear 
Figures_vary;

## Membrane simulation
clear;
VaryNa = 0;			# HH values
Figures_Membrane;

## VaryNa = 1;			# Modified HH values
## Figures_Membrane;
## MV = [MV;mV];
## MV_eq = [MV_eq ;mV_eq];

## figure(40);
## plot(mt,MV);
## Legend = {"Na" sprintf("%g*Na",vary)};
## legend(Legend)
## xlabel("t (msec)");
## ylabel("V(mV)");
## ##axis([0 2 -100 150])
## fig(Name,"VV",2);

## ## Axon simulation
## clear;
## Figures_Axon;




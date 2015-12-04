## Look at varying conentrations
VaryNa = 1;

## Baseline
vary = 1;
Figures_Membrane;
maxE_e_1 =  maxE_e;
maxE_atp_1 = maxE_atp;

## Vary Na_i
VARY = [0.5:0.1:1.5];
MAXE_e = [];
MAXE_atp = [];
mVV = [];
for vary = VARY
    Figures_Membrane;
    MAXE_e = [MAXE_e maxE_e];
    MAXE_atp = [MAXE_atp maxE_atp];
    mVV = [mVV; mV];
endfor

## Two special points
VARYs = [1 1.24];
MAXE_es = [];
MAXE_atps = [];
for vary = VARYs
    Figures_Membrane;
    MAXE_es = [MAXE_es maxE_e];
    MAXE_atps = [MAXE_atps maxE_atp];
endfor

figure(40); 
plot(VARY,MAXE_e/maxE_e_1,"-;Actual;", \
     VARY,MAXE_atp/maxE_e_1,"--;Proxy;",\
     VARYs,MAXE_es/maxE_e_1, "*","linewidth",20,\
     VARYs,MAXE_atps/maxE_e_1, "*","linewidth",20);
grid;
xlabel("Normalised Na^+_i")
ylabel("Normalised Energy")
axis([0.5 1.5 0.7 1.3])
fig("Memb","vary");

Results = [VARY; MAXE_e/maxE_e_1; MAXE_atp/maxE_e_1]

Resultss = [VARYs; MAXE_es/maxE_e_1; MAXE_atps/maxE_e_1]
figure(41); plot(mt,mVV)



## HH & GHK Simulation
names = {"HH" "GHK"};
##names = {"GHK"}
for i_name = 1:length(names)
  name = names{i_name}
  system(sprintf("cp -v HH_%s_numpar.m HH_numpar.m",name));
  clear functions
  for u_HH = [0 20]
    Figures_HH;
  end
endfor

## GHK & MM affinity functions
Ions = {"K" "Na"};
for i_ion = 1:length(Ions)
  Ion = Ions{i_ion};
  Figures_GHK;
endfor


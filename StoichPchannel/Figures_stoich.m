clear;
Setpath;

last = 10;
Set_Membrane;

## Stoichiometry
m = dm2stoich("Membrane","");
n_X = m.n_X;
n_V = m.n_V;

tex = mat2tex(m.N);
printf("N = \\begin{pmatrix}\n %s \n \\end{pmatrix}\n", tex);

## Names
name = "Membrane";
## [input_name,output_name,state_name,nonstate_name] = Membrane_struc;
## species = cellstr((state_name));
## reaction = cellstr((input_name));
species = {"x_m" "n_c" "n_o" "K_i" "K_e" "m_c" "m_o" "h_c" "h_o" "N_i" "N_e" "L_i" "L_e"};
reaction = {"r_n" "r_K" "r_m" "r_h" "r_N" "r_L"};
## Write X 
i_X = [1:m.n_X]';
X_name = species(i_X)

## Write out equations etc.
stoich_info (m.N_f,m.N_r,[],[],name,species,reaction);


for ii = 0:1
  ## Fix
  if ii>0
    i_FIX = i_fix;
  else
    i_FIX = [];
  endif

  n_fix = length(i_FIX);
  Z = zeros(n_fix,n_V);
  N_fix = m.N;
  N_fix(i_FIX,:) = Z;
  N = N_fix;
  m_fix.N = N;
  m_fix = stoich2stoich(m_fix);

  ## Reduced state
  i_x = m_fix.L_xX*i_X;
  x_name = species(i_x)


  ## Conserved moieties
  G = m_fix.G;
  if n_fix==0 # Modify G by hand
     G(4,:) = G(4,:) + G(3,:) - 3*G(1,:);
  endif

## print matrices
  tex = mat2tex(m_fix.N);
  printf("N = \\begin{pmatrix}\n %s \n \\end{pmatrix}\n", tex);
  tex = mat2tex(G);
  printf("G = \\begin{pmatrix}\n %s \n \\end{pmatrix}\n", tex);

  printf("Conserved Moieties\n")
  [n_g m_g] = size(G);
  for i = 1:n_g
    J = find(abs(G(i,:))>1e-6);
    mG = min(abs(G(i,J)));
    printf("%i\t& $",i)
    for j=J
      g_ij = round(G(i,j)/mG);
      if (g_ij>0)
	sig = " + ";
      else
	sig = " - ";
      endif
      if (abs(g_ij)!=1)
	num = sprintf("%i",abs(g_ij));
      else
	num = "";
      endif
      printf("%s%s%s ",sig, num,species{j})
    endfor
    printf("$\\\\\n")
  endfor
endfor

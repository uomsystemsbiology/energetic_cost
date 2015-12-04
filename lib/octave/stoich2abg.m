function stoich2abg (N_f,N_r,Ne_f,Ne_r,name,species,reaction)
	 
  ## usage: stoich2abg (N_f,N_r,Ne_f,Ne_r,name,species,reaction)
  ## 
  ## 

  N = -N_f + N_r;			# Stoichiometric matrix
  NN_f = [N_f;Ne_f];			# Include external species
  NN_r = [N_r;Ne_r];			# Include external species


  ## Sizes
  [n_X,n_V] = size(N_f);	# States
  [nn_X,nn_V] = size(NN_f);	# States + Inputs
  n_U = nn_X - n_X;		# Inputs

  ## Default names
  if nargin<3
     name = "noname";
  endif

  if nargin<4
    for i=1:n_X
	species{i} = sprintf("c_%i",i);
    endfor
  endif

  if nargin<5
    for i=1:n_V
	reaction{i} = sprintf("r_%i",i);
    endfor
  endif

  ## Bond numbers with and without transformer
  n_f1 = sum(sum(NN_f>0));
  ##n_f2 = sum(sum(N_f>1));	

  n_r1 = sum(sum(NN_r>0));
  ##n_r2 = sum(sum(N_r>1));

  # Bond indices
  i_C = [1:n_X];		# C components
  i_SS = [1:n_U] + max(i_C);	# SS components
  i_CSS = [i_C, i_SS];		# C & SS components
  i_f = [1:n_V] + max(i_SS);	# Re components: forward
  i_r = [1:n_V] + max(i_f);	# Re components: reverse
  i_f1 = [1:n_f1] + max(i_r);	# Non-TF forward bonds
  i_r1 = [1:n_r1] + max(i_f1);	# Non-TF reverse bonds
  ## i_f21 = [1:n_f2] + max(i_r1);	# TF forward bonds
  ## i_f22 = [1:n_f2] + max(i_f21);	# TF forward bonds
  ## i_r21 = [1:n_r2] + max(i_f21);	# TF reverse bonds
  ## i_r22 = [1:n_r2] + max(i_f21);	# TF reverse bonds
  n_bonds = max(i_r1);

  ## Label stoichiometric matrices with bond indices.
  ii_f = ii_r = 0;
  II_f = II_r = zeros(nn_X,n_V);
  for i_X = 1:nn_X
    for i_V = 1:n_V
      if NN_f(i_X,i_V)>0
	II_f(i_X,i_V) = i_f1(++ii_f);
      endif
      if NN_r(i_X,i_V)>0
	II_r(i_X,i_V) = i_r1(++ii_r);
      endif
    endfor
  endfor
  
  i_bond = 0;
  i_cmp = 0;

  ## Set default causality
  Causality = zeros(n_bonds,2);	# Initial causality
  ## Effort is imposed by the SS components.
  for i=i_SS
      Causality(i,:) = [-1 -1];
  endfor

  ## Effort is imposed by the C components.
  for i=i_C
      Causality(i,:) = [-1 -1];
  endfor

  ## File
  filename = sprintf("%s_abg.m",name);
  fid = fopen(filename,"w");

  ## Header
  fprintf(fid,"function [%s] = %s_abg\n",name,name);

  fprintf(fid,"subsystemlist = """";\n\n");

  ## C components
  for i_X = 1:n_X
    ##C_name = sprintf("c_%i",i_X);
    C_name = species{i_X};
    fprintf(fid,"\n## C component %s\n", C_name);
    fprintf(fid,"subsystemlist = [subsystemlist; ""%s""];\n",C_name);
    fprintf(fid,"%s.subsystems.%s.type\t = ""C"";\n", name, C_name);
    fprintf(fid,"%s.subsystems.%s.cr\t = ""cp"";\n", name, C_name);
    fprintf(fid,"%s.subsystems.%s.arg\t = ""K_%i,RT"";\n", name, C_name,i_X);
    fprintf(fid,"%s.subsystems.%s.repetitions\t = 1;\n", name, C_name);
    fprintf(fid,"%s.subsystems.%s.status\t = -1;\n", name, C_name);
    fprintf(fid,"%s.subsystems.%s.connections\t = [%i];\n", name, C_name, i_C(i_X));
  endfor

  ## SS components
  for i_U = 1:n_U
    SS_name = species{n_X+i_U};
    fprintf(fid,"\n## SS component %s\n", SS_name);
    fprintf(fid,"subsystemlist = [subsystemlist; ""%s""];\n",SS_name);
    fprintf(fid,"%s.subsystems.%s.type\t = ""SS"";\n", name, SS_name);
    fprintf(fid,"%s.subsystems.%s.cr\t = ""cp"";\n", name, SS_name);
    fprintf(fid,"%s.subsystems.%s.arg\t = ""K_%i,RT"";\n", name, SS_name,i_U);
    fprintf(fid,"%s.subsystems.%s.repetitions\t = 1;\n", name, SS_name);
    fprintf(fid,"%s.subsystems.%s.status\t = -1;\n", name, SS_name);
    fprintf(fid,"%s.subsystems.%s.connections\t = [%i];\n", name, SS_name, i_SS(i_U));
  endfor

  ## Re components
  for i_V = 1:n_V
    Re_name = reaction{i_V};
    fprintf(fid,"\n## Re component %s\n", Re_name);
    fprintf(fid,"subsystemlist = [subsystemlist; ""%s""];\n", Re_name);
    fprintf(fid,"%s.subsystems.%s.type\t = ""Re"";\n", name, Re_name);
    fprintf(fid,"%s.subsystems.%s.cr\t = ""MA"";\n", name, Re_name);
    fprintf(fid,"%s.subsystems.%s.arg\t = ""flow,kappa_%i"";\n", name, Re_name,i_V);
    fprintf(fid,"%s.subsystems.%s.repetitions\t = 1;\n", name, Re_name);
    fprintf(fid,"%s.subsystems.%s.status\t = -1;\n", name, Re_name);
    fprintf(fid,"%s.subsystems.%s.connections\t = [%i -%i];\n", name, Re_name, i_f(i_V), i_r(i_V));
  endfor

  ## 0 components
  for i_X = 1:nn_X
    Z_name = sprintf("z_%s",species{i_X});
    fprintf(fid,"\n## 0 component %s\n", Z_name);
    fprintf(fid,"subsystemlist = [subsystemlist; ""%s""];\n", Z_name);
    fprintf(fid,"%s.subsystems.%s.type\t = ""0"";\n", name, Z_name);
    fprintf(fid,"%s.subsystems.%s.cr\t = """";\n", name, Z_name);
    fprintf(fid,"%s.subsystems.%s.arg\t = """";\n", name, Z_name);
    fprintf(fid,"%s.subsystems.%s.repetitions\t = 1;\n", name, Z_name);
    fprintf(fid,"%s.subsystems.%s.status\t = -1;\n", name, Z_name);    ## i_F = i_f1(find(N_f(i_X,:)==1));
    i_F = setdiff(II_f(i_X,:),[0]);
    i_R = setdiff(II_r(i_X,:),[0]);
    connections = num2str([-i_CSS(i_X) -i_F i_R]);
    fprintf(fid,"%s.subsystems.%s.connections\t = [%s];\n", name, Z_name, connections);
  endfor

  ## 1 components: forward
  for i_V = 1:n_V
    O_name = sprintf("f_%i",i_V);
    fprintf(fid,"\n## 1 component %s\n", O_name);
    fprintf(fid,"subsystemlist = [subsystemlist; ""%s""];\n",  O_name);
    fprintf(fid,"%s.subsystems.%s.type\t = ""1"";\n", name, O_name);
    fprintf(fid,"%s.subsystems.%s.cr\t = """";\n", name, O_name);
    fprintf(fid,"%s.subsystems.%s.arg\t = """";\n", name, O_name);
    fprintf(fid,"%s.subsystems.%s.repetitions\t = 1;\n", name, O_name);
    fprintf(fid,"%s.subsystems.%s.status\t = -1;\n", name, O_name);
    i_F = setdiff(II_f'(i_V,:),[0]);
    connections = num2str([ i_F -i_f(i_V)]);
    fprintf(fid,"%s.subsystems.%s.connections\t = [%s];\n", name, O_name, connections);
  endfor

  ## 1 components: reverse
  for i_V = 1:n_V
    O_name = sprintf("r_%i",i_V);
    fprintf(fid,"subsystemlist = [subsystemlist; ""%s""];\n",  O_name);
    fprintf(fid,"\n## 1 component %s\n", O_name);
    fprintf(fid,"%s.subsystems.%s.type\t = ""1"";\n", name, O_name);
    fprintf(fid,"%s.subsystems.%s.cr\t = """";\n", name, O_name);
    fprintf(fid,"%s.subsystems.%s.arg\t = """";\n", name, O_name);
    fprintf(fid,"%s.subsystems.%s.repetitions\t = 1;\n", name, O_name);
    fprintf(fid,"%s.subsystems.%s.status\t = -1;\n", name, O_name);
    i_R = setdiff(II_r'(i_V,:),[0]);
    connections = num2str([ -i_R i_r(i_V)]);
    fprintf(fid,"%s.subsystems.%s.connections\t = [%s];\n", name, O_name, connections);
  endfor

  ## Subsystem list
  fprintf(fid,"\n## Subsystem list\n");
  fprintf(fid,"%s.subsystemlist = subsystemlist;\n",name);

  ## Create bond list.
  fprintf(fid,"\n## Bonds\n");
  fprintf(fid,"%s.bonds = [\n", name);
  for i_bond=1:n_bonds
    causality = Causality(i_bond,:);
    fprintf(fid,"%s\n", num2str(causality));
  endfor
  fprintf(fid,"];\n");
  fclose(fid);
  
endfunction

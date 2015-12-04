function stoich_info (N_f,N_r,Ne_f,Ne_r,name,species,reaction)
	 
  ## usage: stoich_info (N_f,N_r,Ne_f,Ne_r,name,species,reaction)
  ## 
  ## 
	 
  if nargin<3
    Ne_f = [];
  endif

  if nargin<4
    Ne_r = [];
  endif


  N = -N_f + N_r;			# Stoichiometric matrix
  NN_f = [N_f;Ne_f];			# Include external species
  NN_r = [N_r;Ne_r];			# Include external species

  ## Sizes
  [n_X,n_V] = size(N_f);	# States
  [nn_X,nn_V] = size(NN_f);	# States + Inputs
  n_U = nn_X - n_X;		# Inputs

  ## Default names
  if nargin<5
     name = "noname";
  endif

  if nargin<6
    for i=1:n_X
	species{i} = sprintf("c_%i",i);
    endfor
  endif

  if nargin<7
    for i=1:n_V
	reaction{i} = sprintf("r_%i",i);
    endfor
  endif

  ## Find isolated Cs
  I_C = find(sum(N_f+N_r,2)==0);
  if length(I_C)>0
     printf("The following species are isolated:\n");
     for i=I_C
       printf("\t%s\n", species{i});
     endfor
  endif

  ## Find disconnected reactions
  I_f = find(sum(NN_f,1)==0);
  if length(I_f)>0
     printf("The following reactions have no forward connection:\n");
     for i=I_f
       printf("\t%s\n", reaction{i});
     endfor
  endif

  I_r = find(sum(NN_r,1)==0);
  if length(I_r)>0
     printf("The following reactions have no reverse connection:\n");
     for i=I_r
       printf("\t%s\n", reaction{i});
     endfor
  endif

  ## ODES for species
  filename = sprintf("%s_ode.tex",name);
  fid = fopen(filename,"w");
  for i_X=1:n_X
    fprintf(fid,"\\dot{X}_{%s} &= ",species{i_X})
    for i_V=1:n_V
	if N(i_X,i_V)!=0
	   if N(i_X,i_V)>0
	      pm = "+";
	   else
	     pm = "-";
	   endif
	   if abs(N(i_X,i_V))>1
	      num = num2str(abs(N(i_X,i_V)));
	   else
	       num = "";
	   endif
	  fprintf(fid," %s %sV_{%s} ",pm,num,reaction{i_V});
	endif
    endfor
    #if i_X==n_X
    #  eol = " ";
    #else
      eol = "\\\\";
    #endif
    fprintf(fid,"%s\n",eol);
  endfor
  fclose(fid);

  ## Flow equations
  filename = sprintf("%s_flow.tex",name);
  fid = fopen(filename,"w");
  fr = ["f" "r"];
  for i_V=1:n_V
    fprintf(fid,"V_{%s} \t&= k(X)(",reaction{i_V})
    for i_fr=[1 2];
      for i_X=1:nn_X
	n = eval(sprintf("NN_%s'(i_V,i_X);",fr(i_fr)));
	if n!=0
	  if n>1
	    ex = sprintf("^{%i}",n);
	  else
	    ex = "";
	  endif
	  fprintf(fid,"X_{%s}%s ",species{i_X},ex);
	endif
      endfor
      if i_fr==1
	fprintf(fid," - K_{%s}",reaction{i_V});
      endif
    endfor
    ## if i_V=n_V
    ##   eol = "";
    ## else
      eol = "\\\\";
    ##endif			
    fprintf(fid,")%s\n",eol);
  endfor
  fclose(fid);
  
  ## Reactions
  filename = sprintf("%s_reac.tex",name);
  fid = fopen(filename,"w");
  Eqn_list = reaction2str (1:n_V,species,NN_f,NN_r);
  for i_V = 1:n_V
      fprintf(fid,"\\text{%s}:\t&%s\\\\\n",reaction{i_V}, Eqn_list{i_V});
  endfor
  fclose(fid);

endfunction

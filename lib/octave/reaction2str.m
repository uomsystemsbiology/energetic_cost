function Str = reaction2str (I_v,species,Nf,Nr)
	 
  ## usage: Str = reaction2str (I_v,species,Nf,Nr)
  ## 
  ## 

  [n_n,v_x] = size(Nf);
	 
  i_str = 0;
  for i_v = I_v
    ## Reactants
    NTf = Nf';
    reac = "";
    j_v = find(NTf(i_v,:)>0);
    jj = 0;
    for j=j_v
      if jj++==0
	plus = "";
      else
	plus = "+";
      endif
      n = NTf(i_v,j);
      if n==1
	n_str = "";
      else
	n_str = num2str(n);
      endif
      reac = sprintf("%s%s%s%s", reac, plus, n_str, species{j});
    endfor
    
    ## Products
    NTr = Nr';
    prod = "";
    j_v = find(NTr(i_v,:)>0);
    jj = 0;
    for j=j_v
      if jj++==0
	plus = "";
      else
	plus = "+";
      endif
      n = NTr(i_v,j);
      if n==1
	n_str = "";
      else
	n_str = num2str(n);
      endif
      prod = sprintf("%s%s%s%s", prod, plus, n_str, species{j});
    endfor
    
    str = sprintf("%s = %s", reac,prod);
    Str{++i_str} = str;
  endfor

  ## Return string if only one reaction requested.
  if length(I_v)==1
     Str = str;
  endif

endfunction

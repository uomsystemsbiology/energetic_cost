function [s,p,n_s,n_p] = stoic2reactants (reaction,bg)
	 
  ## usage: [s,p,n_s,n_p] = stoic2reactants (reaction,bg)
  ## 
  ## reaction: reaction name
  ## bg data structure containing N_f and N_r and index.
  ## s cell array of substrate names
  ## p cell array of product names
  ## n_s stoichiometry of substrate
  ## n_p stoichiometry of product

  ## Find reaction index
  ind = eval(sprintf("bg.r_index.%s;", reaction));
  
  ##Find substrate names
  J = find(bg.N_f'(ind,:)>0);
  for i=1:length(J);
    s{i} = bg.species{J(i)};
    n_s(i) = bg.N_f'(ind,J(i));
  endfor

  ##Find product names
  J = find(bg.N_r'(ind,:)>0);
  for i=1:length(J);
    p{i} = bg.species{J(i)};
    n_p(i) = bg.N_r'(ind,J(i));
  endfor
 
endfunction

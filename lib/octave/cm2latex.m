function tex = cm2latex (G,species)
	 
## usage: tex = cm2latex (G,species)
## 
## Tex string coresponding to the conserved moieties

  [n,m] = size(G);
  tex = "";
  for i = 1:n
    for j = 1:m
	if G(i,j)!=0
	  if G(i,j)==1
	    sym = "+";
	  elseif G(i,j)==-1
	    sym = "-";
	  elseif G(i,j)<-1
	    sym = num2str(G(i,j));
	  else
	      sym = sprintf("+%i",G(i,j));
	  endif
	   tex = sprintf("%s %s x_%s", tex, sym, species{j});
	endif
    endfor
    if i<n
       tex = sprintf("%s \\\\\n",tex);
    else
       tex = sprintf("%s \n",tex);
    endif
  endfor
	 
endfunction

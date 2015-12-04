function str = stoich2eqn (M,name,var)
	 
## usage: str = stoich2eqn (M,name,var)
## 
## 

  [n,m] = size(M);
  for i=1:n
      st = "";
      for j=1:m
	if M(i,j)>0
	  st = sprintf("%s + %s_{%s}", st, var, name{j});
	elseif M(i,j)<0
	  st = sprintf("%s -  %s_{%s}", st, var, name{j});
	endif
      endfor
      str{i} = st;
  endfor
	  
endfunction

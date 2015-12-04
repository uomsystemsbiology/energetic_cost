function perm = name2perm (name1,name2,ignore_case)
	 
## usage: perm = name2perm (name1,name2)
## 
## 

	 
  if nargin<3
    ignore_case=1;
  endif

  ## Sanity check
  n = length(name1);
  if (n!=length(name2))
     error("name1 and name2 must be cell arrays with same number of elements");
  endif 

  if ignore_case
     name1 = tolower(name1);
     name2 = tolower(name2);
  endif

  perm = zeros(1,n);
  for i=1:n
      for j=1:n
	  if strcmp(name1{i},name2{j})
	     perm(i) = j;
	  endif
      endfor
  endfor
  
  i_missing = find(perm==0);
  for i = i_missing
      warning(sprintf("Element %i (%s) has no match", i, name1{i})); 
  endfor
endfunction

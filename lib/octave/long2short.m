function short = long2short (long)

  ## usage:  short = long2short (long)
  ##
  ## 

  short = "";
  for i = 1:rows(long)
    l = long(i,:);
    aa = strsplit(strrep(l,"__","%"), "%");
    n_aa = length(aa);
    last = deblank(aa{n_aa});
    if strcmp(last,"u_u")	# Re component
      i_name = n_aa-2;
    else
      i_name = n_aa;
    endif
    if strcmp(aa{i_name},"Re")	# Rea component
       i_name = i_name-1;
    endif
    s = aa{i_name};

    ## Allow for hierarchical names (prefixed by mtt)
    while (length(strfind(s,"mtt")>0)) 
      i_name++;
      s = aa{i_name};
    endwhile

    ## Split at _ to avoid repeated state name
    aa = strsplit(s,"_");
    s = aa{1};

    short = [short;s];
  endfor
  
endfunction

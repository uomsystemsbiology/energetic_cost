function name = mtt_strip_name (fullname,sysname)

  ## usage:  name = strip_name (fullname,sysname)
  ##
  ## Strips the system name from the full name
  
  ## name format is sysname__name
  
  [N,M] = size(fullname);
  
  name = "";
  for i = 1:N
    
    fn = fullname(i,:);

    ## Sanity check
    if !(length(strmatch (sysname,fn))==1)
      error(sprintf("%s does not begin with %s", fn, sysname));
    endif
    
    ## Nasty fudge to fix unnamed systems
    fn = strrep(fn, "mtti", "mtt");
    n_fn = length(fn);
     ## Save up array
    name = [name; fn(length(sysname)+3:n_fn)];
    
  endfor
  
endfunction

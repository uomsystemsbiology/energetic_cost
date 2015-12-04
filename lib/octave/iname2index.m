function ind = iname2index (name,i_name)

  ## usage:  i = iname2index (iname)
  ##
  ## 

  ind = [];
  for i = 1:rows(name)
    ii = strmatch(deblank(name(i,:)),long2short(i_name));
    if length(ii)>0
      ind = [ind;ii];
    end
  end
  

endfunction
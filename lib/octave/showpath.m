function [dr_name neg_name ir_name path] = showpath (w,K,r_name,i_name)

  ## usage: [dr_name path]  = showpath (w,K,r_name,i_name)
  ##
  ## 

  ## Independent reactions
  ir_name = long2short(i_name(find(w!=0),:));

  ## Dependent reactions
  ##  KK = sum(K(:,j),2);		# Consistent weights ????
  KK = K*w;			# Weighted sum
  i_path = find(KK!=0);
  path = KK(i_path);
  dr_name = long2short(r_name(i_path,:));

  ## Reactions with negative sign
  i_neg = find(path<0);
  if length(i_neg)>0
    neg_name = dr_name(i_neg,:);
  else
    neg_name = "";
  end
  
endfunction

function tex = mat2tex(mat)
  ## mat2tex: converts matrix into LaTeX form.
  ## tex = mat2tex(mat)
  

  ###############################################################
  ## Version control history
  ###############################################################
  ## $Id: mat2tex.m,v 1.1 2001/05/10 11:44:40 gawthrop Exp $
  ## $Log: mat2tex.m,v $
  ## Revision 1.1  2001/05/10 11:44:40  gawthrop
  ## Useful conversion functions
  ##
  ## Revision 1.1  1999/03/25 01:33:51  peterg
  ## Initial revision
  ##
  ###############################################################

  [n,m]=size(mat);

  for i=1:n
    for j=1:m
      if mat(i,j) == 0
	printf("0    ");
	else
	  printf("%5.2f", mat(i,j));
	endif
      if j<m
 	printf("&");
      endif
    endfor
      printf("\\\\\n");
  endfor
  
endfunction




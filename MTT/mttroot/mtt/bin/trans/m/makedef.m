function makedef(structure,deffilenum);

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## %% Version control history
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## %% $Id: makedef.m,v 1.9 2001/11/03 07:51:39 geraint Exp $
  ## %% $Log: makedef.m,v $
  ## %% Revision 1.9  2001/11/03 07:51:39  geraint
  ## %% fflush() buffer before returning - ensures def.r is written to.,
  ## %% required for either octave-2.1.35 or linux-2.4.12 (not sure which).
  ## %%
  ## %% Revision 1.8  2000/10/20 13:26:41  peterg
  ## %% Made sure that mttui not declared twice
  ## %%
  ## %% Revision 1.7  2000/10/20 13:16:29  peterg
  ## %% Reformated
  ## %%
  ## %% Revision 1.6  1996/12/07 18:21:57  peterg
  ## %% Now uses fopen + file number
  ## %%
  ## %% Revision 1.5  1996/11/09 21:05:44  peterg
  ## %% Only generates MTTIm when at least 2 states!
  ## %%
  ## %% Revision 1.4  1996/08/30  19:42:36  peter
  ## %% Added newline at end of file.
  ## %%
  ## %% Revision 1.3  1996/08/24 15:06:22  peter
  ## %% Write `END;' at end to please reduce.
  ## %%
  ## %% Revision 1.2  1996/08/18 20:05:20  peter
  ## %% Put unded version control
  ## %%
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



  states = structure(1);
  nonstates=structure(2);
  inputs=structure(3);
  outputs=structure(4);
  zero_outputs = structure(5);
  internal_inputs = structure(6);
  connecting_inputs = structure(7);

  pc = '%';
  ## Declare reduce constants;
  fprintf(deffilenum, 'MTTNx := %1.0f;\n', states);
  fprintf(deffilenum, 'MTTNz := %1.0f;\n', nonstates);
  fprintf(deffilenum, 'MTTNu := %1.0f;\n', inputs);
  fprintf(deffilenum, 'MTTNy := %1.0f;\n', outputs);
  fprintf(deffilenum, 'MTTNyz := %1.0f;\n', zero_outputs);
  fprintf(deffilenum, 'MTTNui := %1.0f;\n', internal_inputs);
  fprintf(deffilenum, 'MTTNuc := %1.0f;\n', connecting_inputs);

  ## Declare reduce matrices
  fprintf(deffilenum, '%s Declare reduce matrices\n', pc);
  if states>0
    fprintf(deffilenum, 'matrix MTTx(%1.0f,1);\n', states);
    fprintf(deffilenum, 'matrix MTTdx(%1.0f,1);\n', states);
  endif

  if nonstates>0
    fprintf(deffilenum, 'matrix MTTz(%1.0f,1);\n', nonstates);
    fprintf(deffilenum, 'matrix MTTdz(%1.0f,1);\n', nonstates);
  endif

  if inputs>0
    fprintf(deffilenum, 'matrix MTTu(%1.0f,1);\n', inputs);
    fprintf(deffilenum, 'matrix MTTdu(%1.0f,1);\n', inputs);
  endif

  if outputs>0
    fprintf(deffilenum, 'matrix MTTy(%1.0f,1);\n', outputs);
  endif

  if zero_outputs>0
    fprintf(deffilenum, 'matrix MTTyz(%1.0f,1);\n', zero_outputs);
    fprintf(deffilenum, 'matrix MTTui(%1.0f,1);\n', zero_outputs);
  elseif internal_inputs>0
    fprintf(deffilenum, 'matrix MTTui(%1.0f,1);\n', inputs);
  endif

  if connecting_inputs>0
    fprintf(deffilenum, 'matrix MTTuc(%1.0f,1);\n', connecting_inputs);
  endif


  ## Make an Nx x Nx unit matrix
  if states>0
    fprintf(deffilenum, 'matrix MTTI(%1.0f,%1.0f);\n', states,states);
    for i = 1:states
      fprintf(deffilenum, 'MTTI(%1.0f,%1.0f) := 1;\n', i, i);
    end
  endif


  ## Make an Nx/2 x Nx/2 unit matrix
  if states>1
    fprintf(deffilenum, 'matrix MTTIm(%1.0f,%1.0f);\n', states/2,states/2);
    for i = 1:states/2
      fprintf(deffilenum, 'MTTIM(%1.0f,%1.0f) := 1;\n', i, i);
    end
  endif


  ## Set the y, yz, du, u, x and dx matrices
  fprintf(deffilenum, '%s Set the y, yz, u and x matrices\n', pc);
  for i=1:outputs
    fprintf(deffilenum, 'MTTy(%1.0f,1) := MTTy%1.0f;\n', i, i);
  endfor

  for i=1:zero_outputs
    fprintf(deffilenum, 'MTTyz(%1.0f,1) := MTTyz%1.0f;\n', i, i);
    fprintf(deffilenum, 'MTTui(%1.0f,1) := MTTui%1.0f;\n', i, i);
  endfor

  for i=1:inputs
    fprintf(deffilenum, 'MTTu(%1.0f,1) := MTTu%1.0f;\n', i, i);
  endfor

  for i=1:inputs
    fprintf(deffilenum, 'MTTdu(%1.0f,1) := MTTdu%1.0f;\n', i, i);
  endfor

  for i=1:states
    fprintf(deffilenum, 'MTTx(%1.0f,1) := MTTx%1.0f;\n', i, i);
  endfor

  for i=1:nonstates
    fprintf(deffilenum, 'MTTdz(%1.0f,1) := MTTdz%1.0f;\n', i, i);
  endfor

  fflush (deffilenum);

endfunction

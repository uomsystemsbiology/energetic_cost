%% Reduce comands to simplify output for system DCexample (DCexample_simp.r)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Version control history
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% $Id: DCexample_simp.r,v 1.1 2000/12/28 17:36:21 peterg Exp $
% %% $Log: DCexample_simp.r,v $
% %% Revision 1.1  2000/12/28 17:36:21  peterg
% %% To RCS
% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% System
%J_m := 1; L_m := 0.1; k_m := 1; b_m := 0.1; r_m := 1;

% Controller - poles at -2 and -10
%alpha_c1 := 12; alpha_c2 := 20;

% Observer -  poles at -10 -50
%alpha_o1 := 60; alpha_o2 := 500;

END;

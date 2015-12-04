%% Reduce comands to simplify output for system Satellite (Satellite_simp.r)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Version control history
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% $Id: Satellite_simp.r,v 1.1 2000/12/28 17:38:40 peterg Exp $
% %% $Log: Satellite_simp.r,v $
% %% Revision 1.1  2000/12/28 17:38:40  peterg
% %% To RCS
% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% System

% Controller - poles -1
alpha_c1 := 2; alpha_c2 := 1;

% Observer -  poles at -1 
alpha_o1 := 2; alpha_o2 := 1;

ON FACTOR;

END;

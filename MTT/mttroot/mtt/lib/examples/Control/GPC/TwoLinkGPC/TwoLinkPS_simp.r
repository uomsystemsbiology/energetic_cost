%% Reduce commands to simplify output for system TwoLinkP (TwoLinkP_simp.r)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Version control history
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% $Id: TwoLinkPS_simp.r,v 1.1 2000/12/28 17:25:33 peterg Exp $
% %% $Log: TwoLinkPS_simp.r,v $
% %% Revision 1.1  2000/12/28 17:25:33  peterg
% %% To RCS
% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LET cos(mttx2) = c_2;
LET sin(mttx2) = s_2;

LET cos(mttx4) = c_4;
LET sin(mttx4) = s_4;

LET cos(mttx2-mttx4) = c;
LET sin(mttx2-mttx4) = s;

END;

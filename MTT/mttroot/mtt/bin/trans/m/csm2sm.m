function [A,B,C,D] = csm2sm(A,B,C,D,E);
% [A,B,C,D] = csm2sm(A,B,C,D,E);
% Converts constrained-state matrices to state matrices.

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Version control history
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% $Id: csm2sm.m,v 1.1 1996/08/19 14:48:48 peter Exp $
% %% $Log: csm2sm.m,v $
% %% Revision 1.1  1996/08/19 14:48:48  peter
% %% Initial revision
% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


A = E\A;
B = E\B;


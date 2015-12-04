function [bonds,status] = ES_cause(bonds);
% Causality for ES component

%SUMMARY ES: converts thermal pseudobond to true bond and vv
%DESCRIPTION Port [e]: temperature (T) / enthalpy flow (E) bond
%DESCRIPTION Port [s]: temperature (T) / entropy flow (S) bond



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Version control history
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% $Id: ES_cause.m,v 1.1 1997/09/04 08:34:53 peterg Exp $
% %% $Log: ES_cause.m,v $
% %% Revision 1.1  1997/09/04 08:34:53  peterg
% %% Initial revision
% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Copyright (c) P.J. Gawthrop, 1997.

status = -1;

% Check that there are exactly two bonds.
if check_bonds(bonds,2,'ES')==0
  return
end

% Same causality as TF
[bonds,status] = TF_cause(bonds);



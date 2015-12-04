function [bonds,status] = RST_cause(bonds);
% Causality for RST component

%SUMMARY RST: thermal resistor with entropy generation
%DESCRIPTION Port[in]: a simple one port resistor but:
%DESCRIPTION flow is multiplied by e_out
%DESCRIPTION Port[out]: thermal power - effort is temperature, flow is entropy flow
%DESCRIPTION Not bicausal (yet)


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Version control history
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% $Id: RST_cause.m,v 1.1 2000/12/28 11:50:45 peterg Exp $
% %% $Log: RST_cause.m,v $
% %% Revision 1.1  2000/12/28 11:50:45  peterg
% %% Initial revision
% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Copyright (c) P.J. Gawthrop, 1997.

status = -1;

% Check that there are exactly two bonds.
if check_bonds(bonds,2,'RST')==0
  return
end

% Same causality as RS
[bonds,status] = RS_cause(bonds);


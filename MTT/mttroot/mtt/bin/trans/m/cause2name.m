function name = cause2name(causality)
% cause2name - converts causality to a string
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%     %%%%% Model Transformation Tools %%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Matlab function  cause2name
% name = cause2name(causality)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Version control history
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% $Id: cause2name.m,v 1.1 2000/12/27 17:00:16 peterg Exp $
% %% $Log: cause2name.m,v $
% %% Revision 1.1  2000/12/27 17:00:16  peterg
% %% Initial revision
% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Copyright (c) P.J. Gawthrop, 1996.

if causality==1
  name = 'effort';
elseif causality ==-1
  name = 'flow';
else
  name = 'state';
end;

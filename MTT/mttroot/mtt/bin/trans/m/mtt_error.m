function mtt_error(message, errorfile);
% mtt_error(message, errorfile);

% function [bonds, status] = abg2cbg(bonds,components,system_name,filename)
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%     %%%%% Model Transformation Tools %%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Version control history
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% $Id: mtt_error.m,v 1.2 2004/09/12 22:27:27 geraint Exp $
% %% $Log: mtt_error.m,v $
% %% Revision 1.2  2004/09/12 22:27:27  geraint
% %% Appended 't' to fopen mode string to open in text mode.
% %%
% %% Revision 1.1  1998/12/03 16:21:27  peterg
% %% Initial revision
% %%
% %% Revision 1.2  1997/02/11 10:06:42  peterg
% %% Removed a debugging line.
% %%
% %% Revision 1.1  1996/08/18  20:06:57  peter
% %% Initial revision
% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set default file if it isn't there already
if nargin<2
  errorfile = fopen('mtt_error.txt','at');
end;

fprintf(errorfile, '*MTT ERROR : %s\n', message);

#error('MTT exiting this transformation on error');


function [bonds,status] = zero_cause(bonds);
% zero_cause - causality for zero junctions
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%     %%%%% Model Transformation Tools %%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Matlab function  zero_cause
% [bonds,status] = zero_cause(bonds);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Version control history
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% $Id: zero_cause.m,v 1.1 2000/12/28 11:50:14 peterg Exp $
% %% $Log: zero_cause.m,v $
% %% Revision 1.1  2000/12/28 11:50:14  peterg
% %% Initial revision
% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Effort
[bonds,e_status] = juncause(bonds,1,1);

%Flow
[bonds,f_status] = juncause(bonds,1,-1);

status = 0;
if (e_status==1)|(f_status==1)       % Over causal
  status = 1;
endif;
if (e_status==-1)|(f_status==-1) % Under causal
  status = -1;
endif;

  

function [bonds,status] = GY_cause(bonds);
% GY_cause - causality for GY component
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%     %%%%% Model Transformation Tools %%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Matlab function  GY_cause
% [bonds,status] = GY_cause(bonds);

%SUMMARY GY: elementary gyrator component
%DESCRIPTION Energy conserving two-port
%DESCRIPTION e_1 = f(f_2); f_1 = f(e_2)



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Version control history
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% $Id: GY_cause.m,v 1.3 1998/05/28 09:59:45 peterg Exp $
% %% $Log: GY_cause.m,v $
% %% Revision 1.3  1998/05/28 09:59:45  peterg
% %% Added a - sign on the first line -- bug fix.
% %%
% %% Revision 1.2  1997/09/12  09:42:12  peterg
% %% Fixed causality bug.
% %%
% %% Revision 1.1  1996/11/01  12:04:25  peterg
% %% Initial revision
% %% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Copyright (c) P.J. Gawthrop, 1996.

%Causality of GY is same as that of a TF but with flipped effort/flow
bonds(2,:) = - bonds(2,2:-1:1);
[bonds,status] = TF_cause(bonds);
bonds(2,:) = - bonds(2,2:-1:1);

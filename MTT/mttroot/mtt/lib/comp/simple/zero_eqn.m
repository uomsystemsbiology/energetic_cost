function structure =  zero_eqn(name, bond_number,bonds,direction,cr,args, ...
    structure,eqnfile);
% zero_eqn - equations for zero junctions
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%     %%%%% Model Transformation Tools %%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Matlab function  zero_eqn
% zero_eqn(name, bond_number,bonds,direction,cr,args, structure,eqnfile);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Version control history
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% $Id: zero_eqn.m,v 1.1 1998/08/25 13:59:48 peterg Exp $
% %% $Log: zero_eqn.m,v $
% %% Revision 1.1  1998/08/25 13:59:48  peterg
% %% Initial revision
% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin<7
  eqnfile = 'stdout';
end;


juneqn(name, bond_number,bonds,direction,cr,args,1,eqnfile)





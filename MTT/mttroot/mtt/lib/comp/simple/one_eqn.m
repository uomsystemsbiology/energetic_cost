function  structure = one_eqn(name,bond_number,bonds,direction,cr,args, ...
    structure,eqnfile);
% one_eqn - equations for a one junction
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%     %%%%% Model Transformation Tools %%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Matlab function  one_eqn
% one_eqn(name,bond_number,bonds,direction,cr,args, ...
%    structure,eqnfile);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Version control history
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% $Id: one_eqn.m,v 1.2 2000/12/28 11:46:33 peterg Exp $
% %% $Log: one_eqn.m,v $
% %% Revision 1.2  2000/12/28 11:46:33  peterg
% %% *** empty log message ***
% %%
% %% Revision 1.1  1996/12/05 10:00:49  peterg
% %% Initial revision
% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp("one-eqn")
name,bond_number,bonds,direction,cr,args
if nargin<7
  eqnfile = 'stdout';
end;

juneqn(name,bond_number,bonds,direction,cr,args,-1,eqnfile)





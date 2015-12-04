function u = mtt_isunique(x);
% if all the elements of x are different, returns 1 else returns 0


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Version control history
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% $Id: mtt_isunique.m,v 1.1 2012/10/15 19:22:44 geraint Exp $
% %% $Log: mtt_isunique.m,v $
% %% Revision 1.1  2012/10/15 19:22:44  geraint
% %% mtt_isunique replaces unique which masks a built-in function.
% %%
% %% Revision 1.1  1997/09/18 19:38:58  peterg
% %% Initial revision
% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = length(x);
u=1

disp("--checking uniqueness --")
x
for i=1:n
  for j = 1:i-1
    if iscell(x)
      if x{i}==x{j}
	u = 0
	break;
      end;
    else
      if x(i)==x(j)
	u = 0
	break;
      end;
    end
  end;
end;


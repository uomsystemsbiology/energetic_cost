function index = name_in_list(name,list)

% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%     %%%%% Model Transformation Tools %%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Matlab function  name_in_list.m
% Finds name in list and returns index -- zero if not found


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Version control history
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% $Id: name_in_list.m,v 1.2 1998/01/23 09:23:00 peterg Exp $
% %% $Log: name_in_list.m,v $
% %% Revision 1.2  1998/01/23 09:23:00  peterg
% %% Fixed typing error in (length(j)==1)
% %%
% %% Revision 1.1  1998/01/23 09:22:05  peterg
% %% Initial revision
% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[N,M] = size(list);
if N<1
  error('list must have at least one element');
end;

index = 0;

for i = 1:N
 j = findstr(name,list(i,:));
 if (length(j)==1)
    if j(1)==1
      if index==0
        index = i;
      else 
        index = [index i];
      end
    end
 end;
end;


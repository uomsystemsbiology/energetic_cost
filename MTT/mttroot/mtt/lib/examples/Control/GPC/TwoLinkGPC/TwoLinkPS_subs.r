%% Reduce  substitution statements for system TwoLinkPSX (TwoLinkPSX_subs.r)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Version control history
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% $Id: TwoLinkPS_subs.r,v 1.1 2000/12/28 17:25:33 peterg Exp $
% %% $Log: TwoLinkPS_subs.r,v $
% %% Revision 1.1  2000/12/28 17:25:33  peterg
% %% To RCS
% %%
% %% Revision 1.2  1998/03/22 20:13:25  peterg
% %% Trig simplification added
% %%
% %% Revision 1.1  1998/03/22 20:12:51  peterg
% %% Initial revision
% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        trig1 := {cos(~x)*cos(~y) => (cos(x+y)+cos(x-y))/2,
                  cos(~x)*sin(~y) => (sin(x+y)-sin(x-y))/2,
                  sin(~x)*sin(~y) => (cos(x-y)-cos(x+y))/2,
                  cos(~x)^2       => (1+cos(2*x))/2,
                  sin(~x)^2       => (1-cos(2*x))/2};
       LET trig1;

END;


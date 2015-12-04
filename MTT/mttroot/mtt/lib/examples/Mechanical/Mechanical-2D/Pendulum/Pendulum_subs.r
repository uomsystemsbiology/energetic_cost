%% Reduce comands to sunstitute in system Pendulum (Pendulum_subs.r)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Version control history
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% $Id: Pendulum_subs.r,v 1.1 2000/12/28 18:01:28 peterg Exp $
% %% $Log: Pendulum_subs.r,v $
% %% Revision 1.1  2000/12/28 18:01:28  peterg
% %% To RCS
% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        trig1 := {cos(~x)*cos(~y) => (cos(x+y)+cos(x-y))/2,
                  cos(~x)*sin(~y) => (sin(x+y)-sin(x-y))/2,
                  sin(~x)*sin(~y) => (cos(x-y)-cos(x+y))/2,
                  cos(~x)^2       => (1+cos(2*x))/2,
                  sin(~x)^2       => (1-cos(2*x))/2};
       LET trig1;

END;

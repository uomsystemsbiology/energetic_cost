%% Parameter file for system twolink (twolink_params.r)
%% This file provides symbolic parameters for simplification

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Version control history
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% $Id: twolink_params.r,v 1.1 1996/12/05 12:18:36 peterg Exp $
% %% $Log: twolink_params.r,v $
% %% Revision 1.1  1996/12/05 12:18:36  peterg
% %% Initial revision
% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Useful trig simplifications


        trig1 := {cos(~x)*cos(~y) => (cos(x+y)+cos(x-y))/2,
                  cos(~x)*sin(~y) => (sin(x+y)-sin(x-y))/2,
                  sin(~x)*sin(~y) => (cos(x-y)-cos(x+y))/2,
                  cos(~x)^2       => (1+cos(2*x))/2,
                  sin(~x)^2       => (1-cos(2*x))/2};
       LET trig1;

% Some simplifications -- see book
j := m*l*l/3;

END;

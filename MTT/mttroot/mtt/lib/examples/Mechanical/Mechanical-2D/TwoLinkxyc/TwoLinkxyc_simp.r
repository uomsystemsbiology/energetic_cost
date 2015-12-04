%% Reduce comands to simplify output for system TwoLinkxyc (TwoLinkxyc_simp.r)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Version control history
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% $Id: TwoLinkxyc_simp.r,v 1.1 2000/12/28 18:02:44 peterg Exp $
% %% $Log: TwoLinkxyc_simp.r,v $
% %% Revision 1.1  2000/12/28 18:02:44  peterg
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

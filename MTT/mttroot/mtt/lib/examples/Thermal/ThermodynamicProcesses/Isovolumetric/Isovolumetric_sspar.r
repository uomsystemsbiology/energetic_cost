% Steady-state parameter file (Isovolumetric_sspar.r)
% Generated by MTT at Wed Mar  4 11:02:40 GMT 1998

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Version control history
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % $Id: Isovolumetric_sspar.r,v 1.1 2000/12/28 18:17:57 peterg Exp $
% % $Log: Isovolumetric_sspar.r,v $
% % Revision 1.1  2000/12/28 18:17:57  peterg
% % To RCS
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Set a pressure of 1 bar
P_0 := 10^5;

% Unit initial volume
V_0 := 1;

% Internal energy
U_0 := P_0*V_0/(gamma_g-1);

% Set initial temperature of 300k
T_0 := 300;

% Deduce the mass of gas
m :=  U_0/(T_0*c_v);

% Entropy
S_0 := U_0/T_0;

% Steady-state states
MTTX1 := 	U_0;         % Isovolumetric_cycle_gas (c)
MTTX2 := 	V_0;         % Isovolumetric_cycle_gas (c)
MTTX3 := 	S_0;         % Isovolumetric_cycle_entropy (3)
MTTX4 := 	V_0;         % Isovolumetric_cycle_volume (3)

% Steady-state inputs
MTTU1 := 	0; % Isovolumetric (Heat)
MTTU2 := 	0; % Isovolumetric (Work)
;;END;

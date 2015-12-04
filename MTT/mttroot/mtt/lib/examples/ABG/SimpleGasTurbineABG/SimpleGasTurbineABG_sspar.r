% Steady-state parameter file (SimpleGasTurbineABG_sspar.r)
% Generated by MTT at Thu Mar 26 16:28:59 GMT 1998

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Version control history
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % $Id: SimpleGasTurbineABG_sspar.r,v 1.1 2000/12/28 16:55:29 peterg Exp $
% % $Log: SimpleGasTurbineABG_sspar.r,v $
% % Revision 1.1  2000/12/28 16:55:29  peterg
% % To RCS
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Find stored mass to give combustion chamber pressure p_3 (at
% temperature t_3
m_c := (p_3*v_c)/(t_3*r);

%Equate pressures
p_4 := p_1;
p_2 := p_3;

%Compute ss temperatures (isentropic)
t_2 := t_1*(p_2/p_1)^alpha;
t_4 := t_3*(p_4/p_3)^alpha;

%Find the steady-state work output
w_0 := c_p*(t_3-t_4) - c_p*(t_2-t_1);

%Compute the corresponding load resistance (to absorb that work)
r_l := w_0/(omega_0)^2;

%Unit mass flow
mdot := 1;

%Corresponding shaft speed
omega_0 := mdot/k;

%Compute shaft inertia to give unit time constant (j_s*r_l)
j_s := r_l;

%Find angular momentum to give shaft speed omega_0
mom_0 :=  omega_0*j_s;


% Steady-state states
MTTX1 := 	mom_0;

% Steady-state inputs - combustion temperature
MTTU1 := 	t_3; % SimpleGasTurbineABG (T3)

;;END;


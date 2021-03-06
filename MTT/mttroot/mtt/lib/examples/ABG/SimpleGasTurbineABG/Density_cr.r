% Constitutive relation file for Density (Density_cr.r)
% Generated by MTT at Wed Mar 11 11:01:28 GMT 1998

OPERATOR Density;

% Ideal gas
FOR ALL R,Temperature,Pressure,Nothing
LET Density(density,ideal_gas,R,effort,3,
	Pressure,effort,1,
	Temperature,effort,2,
	Nothing,flow,3
	) = Pressure/(R*Temperature);

FOR ALL R,Temperature,Pressure,Nothing
LET Density(specific_volume,ideal_gas,R,effort,3,
	Pressure,effort,1,
	Temperature,effort,2,
	Nothing,flow,3
	) = (R*Temperature)/Pressure;

% Incompressible
FOR ALL rho,Temperature,Pressure,Nothing
LET Density(density,incompressible,rho,effort,3,
	Pressure,effort,1,
	Temperature,effort,2,
	Nothing,flow,3
	) = rho;

FOR ALL rho,Temperature,Pressure,Nothing
LET Density(specific_volume,incompressible,rho,effort,3,
	Pressure,effort,1,
	Temperature,effort,2,
	Nothing,flow,3
	) = 1/rho;

END;

function v = rate_dm (vp,vm,v_max,k_v,rho_v)
	 
## usage: v = rate_dm (vp,vm,v_max,k_v,rho_v)
## 
## Computes reaction rate v for the "direct binding modular"
## version
 
v = v_max*(vp-vm)./(k_v + (1-rho_v)*vp + rho_v*vm);
 
	 
endfunction

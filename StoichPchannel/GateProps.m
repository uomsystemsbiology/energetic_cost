function [g,tau] = GateProps (alpha,beta)
	 
## usage: [g,tau] = GateProps (alpha,beta);
## 
## 

## Time constant.
tau = 1./(alpha+beta);

## SS gain
g = alpha./(alpha+beta);

endfunction

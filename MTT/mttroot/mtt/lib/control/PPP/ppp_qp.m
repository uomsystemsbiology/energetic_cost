function [u,U,n_active] = ppp_qp (x,W,J_uu,J_ux,J_uw,Us0,Gamma,gamma,mu,test)

  ## usage:  [u,U] = ppp_qp (x,W,J_uu,J_ux,J_uw,Gamma,gamma)
  ## INPUTS:
  ##      x: system state    
  ##      W: Setpoint vector
  ##      J_uu,J_ux,J_uw: Cost derivatives (see ppp_lin)
  ##      Us0: value of U* at tau=0 (see ppp_lin)
  ##      Gamma, gamma: U constrained by Gamma*U <= gamma 
  ##      mu  Parameter of qp_mu
  ## Outputs:
  ##      u: control signal
  ##      U: control weight vector
  ##
  ## Predictive pole-placement of linear systems using quadratic programming
  ## Use ppp_input_constraint and ppp_output_constraint to generate Gamma and gamma
  ## Use ppp_lin to generate J_uu,J_ux,J_uw
  ## Use ppp_cost to evaluate resultant cost function

  ## Copyright (C) 1999 by Peter J. Gawthrop
  ## 	$Id: ppp_qp.m,v 1.6 2004/10/20 21:58:12 gawthrop Exp $	

  if nargin<9
    mu = 0
  endif

  if nargin<10
    test=0;
  endif
  

  ## Check the sizes
  n_x = length(x);

  [n_U,m_U] = size(J_uu);
  if n_U != m_U
    error("J_uu must be square");
  endif

  [n,m] = size(J_ux);
  if (n != n_U)||(m != n_x)
    error("J_ux should be %ix%i not %ix%i",n_U,n_x,n,m);
  endif


  if length(gamma)>0		# Constraints exist: do the QP algorithm
    ## QP solution for weights U	
##    [U,iterations] = qp_mu(J_uu,(J_ux*x-J_uw*W),Gamma,gamma,mu,[],[],0,test);
    [U,n_active] = qp_hild(J_uu,(J_ux*x - J_uw*W),Gamma,gamma);	# 
##    iterations = 0;

    ##U = qp(J_uu,(J_ux*x - J_uw*W),Gamma,gamma); # QP solution for weights U
    ##U = pd_lcp04(J_uu,(J_ux*x - J_uw*W),Gamma,gamma); # QP solution for weights U
    u = Us0*U;			# Control signal
  else			# Do the unconstrained solution
    ## Compute the open-loop gains
    n_active = 0;
    K_w = J_uu\J_uw;
    K_x = J_uu\J_ux;

    ## Closed-loop control
    U = K_w*W - K_x*x;		# Basis functions weights - U(t)
    u = Us0*U;			# Control u(t)
  endif

endfunction

function [K_x,K_w] = ppp_lin_exact (A,B,C,D,A_u,k_x)

  ## usage:  [K_x,K_w] = ppp_lin_exact (A,B,C,D,A_u,k)
  ##
  ## Computes exact values of the ppp gain matrices:
  ##  where: U = K_w w - K_x x
  ## A,B,C,D State matrices
  ## k_x,k_w       Continuous-time state feedback gains
  ##  where: u = k_w w - k_x x

  [n_x,n_u,n_y] = abcddim(A,B,C,D); # Dimensions

  K_x = ppp_open2closed (A_u,A-B*k_x,k_x);

  ## Steady-state x (and u) for w=1;
  Z = zeros(n_x,1);
  xu_w = [A B;C D]\[Z;1];
  x_w = xu_w(1:n_x);
  u_w = xu_w(n_x+1)
  
  K_w = K_x*x_w;			#
endfunction
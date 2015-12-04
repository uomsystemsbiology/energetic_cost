function [A B C D] = stoich2lin (Nf,Nr,kappa,K,Nef,Ner,Ke,x_ss,x_e)
	 
  ## usage: [A B C D] = stoich2lin (Nf,Nr,kappa,K,Nef,Ner,Ke,x_ss,x_e)
  ## 
  ## 

  warning("stoich2lin is deprecated - please use stoich_lin instead");
	 
  ## Required matrices
  N = -Nf+Nr;	 
  Ne = -Nef+Ner;
  NNf = [Nf;Nef];	 
  NNr = [Nr;Ner];
  KK = [K;Ke]	 		# FIXME
  Kappa = diag(kappa);	 
  ## kNK = Kappa*N'*diag(K);
  ## NkNK = N*kNK;
  ## kNe = Kappa*Ne';
  ## NkNe = N*kNe;

  [nx,nv] = size(Nf);
  [nX,nv] = size(NNf);

  if nargin<9
    x_e = ones(nX-nx,1)
  endif
  
  ## if nargin<8
  ##   ## Steady state
  ##   disp("Computing steady-state")
  ##   [x_ss v_ss] = stoich2ss (Nf,Nr,kappa,K,Nef,Ner,Ke,x_e);
  ##   ## x_ss = -(NkNK\NkNe*Ke);
  ##   ## v_ss = -kNK*x_ss - kNe*Ke	 
  ## endif

  ## Steady-state flows
  X_ss = [x_ss;x_e];
  vplus  = Kappa*exp(NNf'*log(KK.*X_ss));
  vminus = Kappa*exp(NNr'*log(KK.*X_ss));
  v_ss = vplus - vminus;

  ## Linearise about steady state
#  C0 = -kNK;
  C = (diag(vplus)*Nf' - diag(vminus)*Nr')/diag(x_ss);
 # CHECKC = norm(C-C0)
  D = (Kappa\diag(v_ss));
  A = N*C;
  B = N*D;
  
endfunction

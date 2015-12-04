function [subport,n]=split_port(port_name);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Version control history
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % $Id: split_port.m,v 1.2 2012/10/15 19:24:19 geraint Exp $
% % $Log: split_port.m,v $
% % Revision 1.2  2012/10/15 19:24:19  geraint
% % split is obsolete in Octave, replaced with strsplit
% %
% % Revision 1.1  1998/04/16 14:08:00  peterg
% % Initial revision
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



  subport = "";
  raw_subport = char(strsplit(port_name, ',')); # Find the components of the vector
				# port 
  [n,m] = size(raw_subport);	# Number of ports

  for i = 1:n
    s = deblank(raw_subport(i,:));
    l = length(s);

    if s(1)=="["		# Remove leading [ (if any)
      l = l-1;
      s = substr(s,2,l);
    end

    if s(l)=="]"		# Remove trailing ] (if any)
      l = l-1;
      s = substr(s,1,l);
    end

    s = ["[",s,"]"];		# Surround by []
    
    subport = [subport ;s];	# Return to the array
  end;

endfunction;






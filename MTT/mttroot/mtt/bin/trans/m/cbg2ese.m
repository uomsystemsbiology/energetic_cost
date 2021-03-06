function structure = cbg2ese(system_name, system_type, system_cr, ...
			     system_args, full_name, full_name_repetition, ...
			     repetition,...
			     structure, structure_file,infofilenum)
  ## Set up globals to count the component inputs and outputs. This relies on
  ## the named SS (the ports) being in the correct order. Using globals here
  ## avoids changing the common argument list for all _eqn files for something
  ## which is only used for named SS components.
  global local_u_index
  global local_y_index
  global at_top_level
  ## 
  ##     ###################################### 
  ##     ##### Model Transformation Tools #####
  ##     ######################################
  ## 
  ## Matlab function  cbg2ese.m
  ## Acausal bond graph to causal bond graph: mfile format
  ## Structure matrix [states,nonstates,inputs,outputs,zero_outputs]
  
  ## ###############################################################
  ## ## Version control history
  ## ###############################################################
  ## ## $Id: cbg2ese.m,v 1.54 2009/11/02 16:54:03 geraint Exp $
  ## ## $Log: cbg2ese.m,v $
  ## ## Revision 1.54  2009/11/02 16:54:03  geraint
  ## ## Replaced deprecated functions from Octave 2.1 for Octave 3.0: is_struct -> isstruct, struct_contains -> isfield, struct_elements -> fieldnames, is_complex -> iscomplex, setstr -> char
  ## ##
  ## ## Revision 1.53  2005/03/21 11:09:47  gawthrop
  ## ## Now handles bicausal SS component -
  ## ##   ie source-source or sensor-sensor
  ## ##
  ## ## Revision 1.52  2004/09/12 22:27:27  geraint
  ## ## Appended 't' to fopen mode string to open in text mode.
  ## ##
  ## ## Revision 1.51  2003/05/16 11:16:28  gawthrop
  ## ## Fixed bug with multiports
  ## ##
  ## ## Revision 1.50  2003/05/08 18:47:50  gawthrop
  ## ## Fixed __ bug when using * representations
  ## ##
  ## ## Revision 1.49  2003/03/13 15:19:04  gawthrop
  ## ## Now uses __ to delimit subsystems in names.
  ## ##
  ## ## Revision 1.48  2003/03/13 15:10:26  gawthrop
  ## ## Removed redundant final column
  ## ##
  ## ## Revision 1.47  2003/02/28 09:12:17  gawthrop
  ## ## Two more columns in _stuc.txt: causality and subsystem name
  ## ##
  ## ## Revision 1.46  2002/08/20 15:51:17  gawthrop
  ## ## Update to work with ident DIY rep
  ## ##
  ## ## Revision 1.45  2002/05/22 09:15:03  gawthrop
  ## ## Non-repetitive components no longer use _1 in names
  ## ##
  ## ## Revision 1.44  2001/11/11 18:12:30  geraint
  ## ## Moved fflush(structure_file) out of loop.
  ## ##
  ## ## Revision 1.43  2001/11/11 08:32:00  geraint
  ## ## fflush (structure_file).
  ## ##
  ## ## Revision 1.42  2001/04/23 16:23:30  gawthrop
  ## ## Now stips ; from bottlom level argument list - allows aliasing of
  ## ## parts of a,b,c (eg a,b by using a,b;c
  ## ##
  ## ## Revision 1.41  2001/04/15 21:15:41  geraint
  ## ## Added interface definition rep: _ICD.(txt|c|cc|m).
  ## ##
  ## ## Revision 1.40  2001/02/05 01:50:29  geraint
  ## ## No unit type comparison at ports if either is "none".
  ## ##
  ## ## Revision 1.40  2000/12/16 08:10:55  geraint
  ## ## No unit type comparison at ports if either is "none".
  ## ##
  ## ## Revision 1.39  2000/11/16 12:54:14  peterg
  ## ## Added checking of unit consistency at ports
  ## ##
  ## ## Revision 1.38  2000/11/16 10:00:57  peterg
  ## ## *** empty log message ***
  ## ##
  ## ## Revision 1.37  2000/11/12 16:45:57  peterg
  ## ## Close ese file before recursive call of cbg2ese -- reopen when
  ## ## finished.
  ## ## THis prevents a new file being opened for each subsystem which fails
  ## ## when > 1K files opened
  ## ##
  ## ## Revision 1.36  2000/10/13 10:54:47  peterg
  ## ## Now writes out a unique name for each state etc
  ## ##
  ## ## Revision 1.35  2000/10/12 19:27:47  peterg
  ## ## Now writes the aliased args
  ## ##
  ## ## Revision 1.34  2000/09/01 08:42:44  peterg
  ## ## Cahged somes ends to end for etc for clarity
  ## ##
  ## ## Revision 1.33  2000/09/01 08:05:32  peterg
  ## ## Reformatted with out changing function
  ## ##
  ## ## Revision 1.32  1998/11/16 13:01:19  peterg
  ## ## Fixed problem with repetitions>1 due to new data structures -- now
  ## ## computes initial next_bond correctly
  ## ##
  ## ## Revision 1.31  1998/09/24 12:57:44  peterg
  ## ## Now ignores aliasing if no arguments given.
  ## ##
  ## ## Revision 1.30  1998/09/02 11:14:23  peterg
  ## ## Revised to use ordered lists of subsystems and ports
  ## ##
  ## ## Revision 1.29  1998/08/25 09:22:34  peterg
  ## ## Correctely recognises port SSs its now easy -- they are in there own
  ## ## field
  ## ##
  ## ## Revision 1.28  1998/08/25 08:31:42  peterg
  ## ## Fixed bug - didn't find the ports - use deblank
  ## ##
  ## ## Revision 1.27  1998/08/25 07:16:49  peterg
  ## ## Modified to data struture representation
  ## ##
  ## ## Revision 1.26  1998/08/24 14:53:55  peterg
  ## ## Uses new _cbg structure.
  ## ##
  ## ## Revision 1.25  1998/07/28 19:05:12  peterg
  ## ## Sttill has vector SS port bug?
  ## ##
  ## ## Revision 1.24  1998/07/27 10:26:02  peterg
  ## ## No change - but fixed bug in alias_args
  ## ##
  ## ## Revision 1.23  1998/07/08 12:33:51  peterg
  ## ## Reinstated the infofilenum parameter.
  ## ##
  ## ## Revision 1.22  1998/07/04 07:10:27  peterg
  ## ## Don't evaluate alias if no constitutive relationship or args and write
  ## ## message.
  ## ##
  ## ## Revision 1.21  1998/07/03 18:58:58  peterg
  ## ## Put arg alias stuff within function alias_args
  ## ## Called recursively to handle arithmetic expressions
  ## ##
  ## ## Revision 1.20  1998/07/03 14:39:09  peterg
  ## ## Added info messages a bit busy now!
  ## ##
  ## ## Revision 1.19  1998/04/12 11:58:19  peterg
  ## ## Rename port components by changing name_r to [name_r
  ## ##
  ## ## Revision 1.18  1998/04/11 18:59:16  peterg
  ## ## at_top_level now global - passed to SS components
  ## ##
  ## ## Revision 1.17  1998/04/04 10:47:31  peterg
  ## ## Uses (coerced) components from _cbg file - _abg not now used here.
  ## ##
  ## ## Revision 1.16  1998/03/06 09:38:58  peterg
  ## ## Now writes correct structure file for multiport C & I
  ## ##
  ## ## Revision 1.15  1997/12/16 18:24:33  peterg
  ## ## Added unknown_input to structure list
  ## ##
  ## ## Revision 1.14  1997/08/26 07:51:10  peterg
  ## ## Now counts the local input and outputs by order of appearence rather
  ## ## than by port number - it therfore handles ports with bicausality correctely.
  ## ##
  ## ## Revision 1.13  1997/05/19 16:45:56  peterg
  ## ## Fixed ISW bug -- deleted spurious ISW_eqn.m file
  ## ##
  ## ## Revision 1.12  1997/04/15  09:17:26  peterg
  ## ## Added the structure file - contains details of states etc.
  ## ## 
  ## ## Revision 1.11  1996/12/07 18:20:11  peterg
  ## ## Replaces null argument by a default and tells user.
  ## ##
  ## ## Revision 1.10  1996/12/07 17:37:07  peterg
  ## ## Now handles numbered SS ports appearing at top level.
  ## ##
  ## ## Revision 1.9  1996/12/04 21:49:47  peterg
  ## ## Compares full-name with empty string (instead of testing for zero
  ## ## length)
  ## ##
  ## ## Revision 1.8  1996/08/30  16:36:08  peter
  ## ## More info written to ese files.
  ## ##
  ## ## Revision 1.7  1996/08/30 11:23:13  peter
  ## ## Argument substitution implemented.
  ## ##
  ## ## Revision 1.6  1996/08/27  08:04:52  peterg
  ## ## Handles complex components and repetative components.
  ## ##
  ## ## Revision 1.5  1996/08/24  15:02:23  peter
  ## ## Writes `END;' to keep reduce happy.
  ## ##
  ## ## Revision 1.4  1996/08/19 09:03:41  peter
  ## ## Handles repeating components.
  ## ##
  ## ## Revision 1.3  1996/08/18 20:08:02  peter
  ## ## Included additional structure: structure(5) = zero_outputs.
  ## ##
  ## ## Revision 1.2  1996/08/08 18:08:11  peter
  ## ## Sorted out file naming sceme
  ## ##
  ## ## Revision 1.1  1996/08/08 15:53:23  peter
  ## ## Initial revision
  ## ##
  ## #############################################################
  
  ## disp("cbg2ese");
  ## system_name, system_type, full_name, repetition
  
  pc = "%";
  sub_delim = "__";		# Subsystem delimiter
  

  unit_error = "Component %s connects inconsistent ports with units %s and %s";
  unit_info = "Component %s connects ports with units %s and %s";  

  ## Set up the names corresponding to the structure matrix.
  structure_name = [
		    "state        ",
		    "nonstate     ",
		    "input        ",
		    "output       ",
		    "zero_output  ",
		    "unknown_input"];
  
  
  ## Are we at the top level of the heirarchy?
  at_top_level = strcmp(full_name, "");
  
  ## Create the (full) system name
  if at_top_level
    full_name = system_name;
    full_name_repetition = system_name;
    system_type = system_name;
  else
    full_name = [full_name, sub_delim, system_name];

    if (repetition>1)
      full_name_repetition = [full_name_repetition, \
			      sub_delim, system_name, sub_delim, \
			      num2str(repetition)];
    else
      full_name_repetition = [full_name_repetition, \
			      sub_delim, system_name];
    endif
    
  end;
  
  
  
  cbg_name = [full_name, "_cbg"];
  if exist(cbg_name)~=2		# Return if cbg file doesn't exist
    disp([cbg_name, " does not exist"]);
    return
  end;
  
  ## Setup files
  ese_name = [full_name_repetition, "_ese.r"];
  ese_file = fopen(ese_name, "wt") # open file (first time)
  icd_file = fopen([full_name_repetition,"_icd.txt2"],"wt")

  fprintf(ese_file, "\n%s%s Equation file for system %s (file %s)\n", ...
	  pc, pc, full_name_repetition, ese_name);
  fprintf(ese_file, "%s%s Generated by MTT\n\n", pc, pc);
  
  ## Evaluate the system function to get the bonds
  eval(["CBG = ", cbg_name, ";"]);
  ##eval(["[bonds,status,system_type,components] = ", cbg_name, ";"]);
  ##  abg_name = [system_type, "_abg"];
  ##  cmp_name = [system_type, "_cmp"];
  ##  alias_name = [system_type, "_alias"];
  
  ## No longer needed - cbg now has components:  
  ##  eval(["[junk,components]=", abg_name, ";"]);
  
  ## Find number of bonds
  [n_bonds,columns] = size(CBG.bonds);
  if (columns ~= 2)&(n_bonds>0)
    error("Incorrect bonds matrix: must have 2 columns");
  endif;

  ## Set up initial bond units
  for i=1:n_bonds
    bond_effort_unit(i,:)="null";
    bond_flow_unit(i,:)="null";
  endfor

  ##  ## Find number of components
  ##  [n_components,columns] = size(components);
  ##  n_components = n_components
  
  ## Set up the first dummy bond number - needed for repetative components
  ##  next_bond = max(max(abs(components)))+1;
  next_bond = n_bonds+1;
  
  ## Set up the counters for the labelled SS. These are, by definition,
  ## encountered first and so the counters will not be messed up by subsystems.
  local_u_index = 0;
  local_y_index = 0;
  
  if (length(system_args)==0)
    mtt_info(sprintf("No arguments given so no argument aliasing done for system %s(%s)",\
		     system_name,system_type), infofilenum);
    AliasingArguments=0;
  else
    AliasingArguments=1;
  endif;
  
  if (length(system_cr)==0)
    mtt_info(sprintf("No cr given so no cr aliasing done for system %s(%s)",\
		     system_name,system_type), infofilenum);
    AliasingCRs=0;
  else
    AliasingCRs=1;  endif;
  
  
  fields=["ports";"subsystems"]; # Do for both ports and subsystems -
  ## ports first
  lists=["portlist";"subsystemlist"];
  for i=1:2
    field=deblank(fields(i,:));
    list=deblank(lists(i,:));
    if isfield(CBG,list);
      eval(["namelist=CBG.",list,";"]); # List of ports/subsystems
      [N,M]=size(namelist);	# Number of ports/subsystems
      for j=1:N
      	comp_name = deblank(namelist(j,:)); # Name of this component
	eval(["subsystem=CBG.",field,".",comp_name,";"]); # Pluck out the details
    	comp = subsystem.connections; # Connections
    	bond_list = abs(comp);
    	direction = sign(comp)'*[1 1];
				# Convert from arrow orientated to component orientated causality
    	comp_bonds = CBG.bonds(bond_list,:).*direction;
	
    	## disp(["---- ", field, " ---"]);    
    	
	if AliasingArguments	# Alias the args list if appropriate
    	  message = sprintf("\tfor component  %s (%s) within %s",\
			    comp_name,subsystem.type,full_name);    
    	  if isfield(CBG,"alias")
	    subsystem.arg = alias_args(subsystem.arg,CBG.alias,";",message,infofilenum,full_name);
    	  endif;
	endif;
	
	if AliasingCRs	# Alias the CR list if appropriate
    	  message = sprintf("\tfor component  %s (%s) within %s",\
			    comp_name,subsystem.type,full_name);    
    	  if isfield(CBG,"alias")
	    subsystem.cr = alias_args(subsystem.cr,CBG.alias,";",message,infofilenum,full_name);
    	  endif;
	endif;

	## Substitute positional ($1 etc) arguments
    	subsystem.cr = subs_arg(subsystem.cr,system_cr, ...
				"lin",full_name,subsystem.type,comp_name,infofilenum);
    	subsystem.arg = subs_arg(subsystem.arg,system_args, ...
				 "1",full_name,subsystem.type,comp_name,infofilenum);
    	
	## change name of 0 and 1 components -- matlab doesn't like numbers here
    	if strcmp(subsystem.type,"0")
	  subsystem.type = "zero";
    	endif;
    	if strcmp(subsystem.type,"1")
	  subsystem.type = "one";
    	endif;
    	
    	ports = length(bond_list);
    	
    	if subsystem.repetitions>1
	  port_pairs = ports/2;
	  if round(port_pairs)~=port_pairs;
	    mtt_info(["Repeated component ", comp_name, ...
		      " has an odd number of ports - ignoring repetitions"], infofilenum);
	    subsystem.repetitions = 1;
	  endif;
    	endif;
    	
    	if subsystem.repetitions>1
	  odd_bonds = bond_list(1:2:ports-1);
	  even_bonds = bond_list(2:2:ports);
    	endif;
    	
    	for k = 1:subsystem.repetitions
	  
	  if subsystem.repetitions>1
	    
	    if k==1
	      bond_list(1:2:ports-1) = odd_bonds;
	    else
	      bond_list(1:2:ports-1) = bond_list(2:2:ports);
	    endif;
	    
	    if k==subsystem.repetitions
	      bond_list(2:2:ports) = even_bonds;
	    else
	      new_bonds = [next_bond:next_bond+port_pairs-1];
	      next_bond = next_bond+port_pairs;
	      bond_list(2:2:ports) = new_bonds;
	    endif;
	    
	  endif;
	  
	  ## Invoke the appropriate equation-generating procedure
	  name_r = full_name_repetition;
	  eqn_name = [subsystem.type, "_eqn"];
	  
	  if exist(eqn_name)~=2 ## Try a compound component
            fclose(ese_file);	# Close but reopen later

	    disp("---PUSH---"); bond_list, comp_name, subsystem.type
	    structure = cbg2ese(comp_name, subsystem.type, subsystem.cr, subsystem.arg, ...
				full_name, full_name_repetition, ...
				k, structure,  structure_file, infofilenum);
	    
	    disp("---POP---");
	    ese_file = fopen(ese_name, "at") # open file (again)

	    eval(["subABG = ",subsystem.type , "_abg;"]); # Get the information

	    ## Link up the bonds for this compound component
	    fprintf(ese_file, ...
		    "\n\t%s Equations linking up subsystem %s (%s)\n\n", ...
		    pc, comp_name, subsystem.type);
	    
	    if (k>1)
	      name_comp_name = sprintf("%s%s%s%s%d", ...
				       full_name_repetition, sub_delim, \
				       comp_name, sub_delim, k);
	    else
	      name_comp_name = sprintf("%s%s%s", ...
				       full_name_repetition, sub_delim, \
				       comp_name);
	    endif
	    
	    
	    printf("\n\t%s Equations linking up subsystem %s (%s)\n\n",\
		   pc, comp_name, subsystem.type);
	    
	    u_index = 0; y_index = 0; ## Count the inputs and outputs
	    for port_number=1:length(bond_list)
	      repetition,port_number
              port_bond_number = bond_list(port_number)
# 	      this_bond_effort_unit = \
# 		  deblank(bond_effort_unit(port_bond_number,:))
# 	      this_bond_flow_unit = \
# 		  deblank(bond_flow_unit(port_bond_number,:));

# 	      ## Extract the unit/domain stuff
#               this_port_name = subABG.portlist(port_number,:);
              
#               eval(sprintf("this_port = subABG.ports.%s;", \
# 			   this_port_name));
# 	      if isfield(this_port,"units")
#                 eval(["effort_unit = \
# 		    subABG.ports.",this_port_name,".units.effort;"]);
#                 eval(["flow_unit = \
# 		    subABG.ports.",this_port_name,".units.flow;"]);
# 	      else
# 		effort_unit = "none";
# 		flow_unit = "none";
# 	      endif

# 	      ## and check consistency
#               ## Efforts
# 	      if strcmp(this_bond_effort_unit,"null") # set
# 		bond_effort_unit = \
# 		    [bond_effort_unit(1:port_bond_number-1,:)
# 		     effort_unit
# 		     bond_effort_unit(port_bond_number+1:n_bonds,:)
# 		     ]
# 	      elseif (!strcmp(this_bond_effort_unit,"none") && !strcmp(effort_unit,"none")) # check
# 		mtt_info(sprintf(unit_info,full_name, effort_unit, \
# 				 this_bond_effort_unit), infofilenum);
# 		if !strcmp(this_bond_effort_unit,effort_unit)
# 		  error_string = sprintf(unit_error, full_name,\
# 					 effort_unit, \
# 					 this_bond_effort_unit);
# 		  mtt_error(error_string);
# 		endif
# 	      endif
# 	      ## Flows
# 	      if strcmp(this_bond_flow_unit,"null") # set
# 		bond_flow_unit = \
# 		    [bond_flow_unit(1:port_bond_number-1,:)
# 		     flow_unit
# 		     bond_flow_unit(port_bond_number+1:n_bonds,:)
# 		     ]
# 	      elseif (!strcmp(this_bond_flow_unit,"none") && !strcmp(flow_unit,"none")) # check
# 		mtt_info(sprintf(unit_info,full_name, flow_unit, \
# 				 this_bond_flow_unit), infofilenum);
# 		if !strcmp(this_bond_flow_unit,flow_unit)
# 		  error_string = sprintf(unit_error, full_name,\
# 					 flow_unit, \
# 					 this_bond_flow_unit);
# 		  mtt_error(error_string);
# 		endif
# 	      endif
	      

	      ## Effort
	      if comp_bonds(port_number,1)==1 # Source
	     	u_index = u_index + 1;
	     	fprintf(ese_file, "%s_MTTu%d := %s;\n", ...
		 	name_comp_name, u_index, varname(name_r, ...
							 bond_list(port_number),1));
	      else # Sensor
	     	y_index = y_index + 1;
	     	fprintf(ese_file, "%s := %s_MTTy%d;\n", ...
		 	varname(name_r, ...
				bond_list(port_number),1), name_comp_name, y_index);
	      end;
	      ## Flow
	      if comp_bonds(port_number,2)==-1 # Source
	     	u_index = u_index + 1;
	     	fprintf(ese_file, "%s_MTTu%d := %s;\n", ...
		 	name_comp_name, u_index, varname(name_r, ...
							 bond_list(port_number),-1));
	      else # Sensor
	     	y_index = y_index + 1;
	     	fprintf(ese_file, "%s := %s_MTTy%d;\n", ...
		 	varname(name_r, ...
				bond_list(port_number),-1), name_comp_name, y_index);
	      end;	
	    end;
	    
	  else # its a simple component
	    fprintf(ese_file, "\n\t%s Equations for component %s (%s), repetition %d\n\n", ...
		    pc, comp_name, subsystem.type,k);
	    
	    
	    ##				# Is it a named port? (Name begins with [)
	    ##            Named_Port = (comp_name(1)=="[");
	    
	    if strcmp(field,"ports") #Add [ to start of name
	      name_r = ["[" name_r];
	    endif;
	    
	    ## Save the current structure
	    old_structure = structure;
	    
	    ## Generate the simple component equations
	    ## .. firstly replacing ; by , in argument list
	    subsystem_arg = strrep(subsystem.arg,";",",");

	    eval(["structure = ", ...
		  eqn_name, ...
		  "(name_r,bond_list,comp_bonds, ...
		    direction,subsystem.cr,subsystem_arg,structure,ese_file);" ]);
	    
	    ## If structure has changed, write info to structure file.
	    structure_change = structure-old_structure;
	    
	    ## The following line is to avoid probs with multiport C or I
	    ## it needs changing to handle multi ports correctly
	    structure_changes = sum(structure_change);
	    
	    structure_change = structure_change>zeros(size(structure_change));
 	    if structure_changes>0
 	      which_indices = getindex(structure_change,1);
 	      which_indices = which_indices(:,2)';
 	      for which_index=which_indices
 	    	value = structure(which_index);
	    	value_change=value-old_structure(which_index);
 	    	for k=1:value_change
		  if strcmp(subsystem.type,"SS") # One port, may be bicausal
		    cause_name = cause2name(-comp_bonds(1,k));
		  else		# Maybe multiport - but unicausal
		    cause_name = cause2name(-comp_bonds(k,1));
		  endif
 		  fprintf(structure_file, ...
 			  "%s\t%i\t%s\t%s%s%s\t%i\t%s\n", ...
 			  structure_name(which_index,:), value-k+1, ...
 			  comp_name, full_name_repetition, sub_delim, comp_name, \
 			  repetition, cause_name); 
 	    	endfor;
 	      endfor;
 	    endif;
	  endif
	endfor
	fflush (structure_file);

	## component interface definition
	if isfield(CBG,"icd")
	  if isfield(CBG.icd,comp_name)
	    if AliasingArguments
	      subsystem.icd = alias_args(eval(["CBG.icd.",comp_name]),CBG.alias,";",message,infofilenum,full_name);
	    endif
	  endif
	  if (isfield(subsystem,"icd"))
	    subsystem.icd = subs_arg(subsystem.icd,system_args, ...
				     "null",full_name,subsystem.type,comp_name,infofilenum);
	    
	    fprintf(icd_file,"%s_%s\t%s\t",full_name_repetition,comp_name,subsystem.icd);
	    if (comp_bonds(1) == 1)
	      fprintf(icd_file,"output,");
	    elseif (comp_bonds(1) == -1)
	      fprintf(icd_file,"input,");
	    endif
	    if (comp_bonds(2) == 1)
	      fprintf(icd_file,"input\n");
	    elseif (comp_bonds(2) == -1)
	      fprintf(icd_file,"output\n");
	    endif
	  endif
	endif			# End of component interface definition

      endfor			# [subsystem,comp_name] = CBG_field
    endif			# isfield(CBG,field)
  endfor			# i=1:2

  fclose(icd_file);
  fclose(ese_file);		# Close
  
endfunction




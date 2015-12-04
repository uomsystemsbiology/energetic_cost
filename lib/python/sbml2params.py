#!/usr/bin/env python

## Required libraries
import sys
import os.path
from libsbml import *

## Main prog
def main (args):
  """Usage: printUnits filename
  """

  if (len(args) != 2):
      print("Usage: printUnits filename");
      return 1;
  
  filename = args[1];
  document = readSBML(filename);
  
  if (document.getNumErrors() > 0):
      print("Encountered the following SBML errors:");
      document.printErrors();
      return 1;
  
  model = document.getModel();

  # Global parameters
  print("## Global parameters");

  for i in range(0,model.getNumParameters()):
      p = model.getParameter(i);
      Name = p.getName();
      names = Name.split("_");
      if (len(names)>1):
          print("{0}.{1} = {2}; # {3}".format(names[0],names[1],p.getValue(),p.getUnits()));

  # Local parameters
  print("## Local parameters");
  reactions = model.getListOfReactions()
  i = 0;
  for reaction in reactions:
      i = i+1;
      rName = reaction.getName();
      print ("## Reaction {0}".format(rName));
      print("reaction{{{0}}} = \"{1}\";".format(str(i),rName));
      paramL = reaction.getKineticLaw().getListOfParameters();
      for p in paramL:
          print("{0}.{1} = {2}; #{3}".format(reaction.getName(), p.getName(), p.getValue(), p.getUnits())); 

  return 0;
  
if __name__ == '__main__':
  main(sys.argv)  

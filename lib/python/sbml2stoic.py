#! /usr/bin/python

## This converts an sbml model to the stoichiometric matrix N and writes to standard output.
## And some other useful info as well.
## Output in Octave format.
## Code hacked from http://www.snip2code.com/Snippet/114023/Parse-SBML-stoichiometry-matrix

import libsbml
from libsbml import *
import argparse
     
def _parser():
    
    parser = argparse.ArgumentParser(description="Parse stoichiometry matrix of SBML file")
    parser.add_argument('file', metavar="filename", type=argparse.FileType('r'),
                        help="Filename of SBML file to parse")
    
    return parser

def _main():
    parser = _parser()
    args = parser.parse_args()
    
    file_ = args.file
    species, reactions, stoichiometry_matrix = parse_file(file_)
    
def parse_file(open_file_):
    reader = libsbml.SBMLReader()
    document = reader.readSBML(open_file_.name)
        
    sbml_model = document.getModel()

    n_species = sbml_model.getNumSpecies()
    n_reactions = sbml_model.getNumReactions()
    print "sbml.N = zeros(", n_species, ",", n_reactions, ");"
    print "sbml.Nf = zeros(", n_species, ",", n_reactions, ");"
    print "sbml.Nr = zeros(", n_species, ",", n_reactions, ");"

    print "sbml.N_float = zeros(", n_species, ",", n_reactions, ");"
    print "sbml.Nf_float = zeros(", n_species, ",", n_reactions, ");"
    print "sbml.Nr_float = zeros(", n_species, ",", n_reactions, ");"

    ## Stoichiometric matrix
    species = [s.getName() for s in sbml_model.getListOfSpecies()]
    reactions = [r.getName() for r in sbml_model.getListOfReactions()]
    stoichiometry_matrix = {}
        
    for reaction_index, reaction in enumerate(sbml_model.getListOfReactions()):
        reactants = {r.getSpecies(): r.getStoichiometry() for r in reaction.getListOfReactants()}
        products = {p.getSpecies(): p.getStoichiometry() for p in reaction.getListOfProducts()}
 
        for species_index, species_node in enumerate(sbml_model.getListOfSpecies()):
            species_id = species_node.getId()
            stoic_f_float = reactants.get(species_id, 0);
            stoic_r_float = products.get(species_id, 0);
            stoic_f = int(reactants.get(species_id, 0));
            stoic_r = int(products.get(species_id, 0));
            ##net_stoichiometry = -int(reactants.get(species_id, 0)) + int(products.get(species_id, 0))
            net_stoichiometry = stoic_r - stoic_f;
            net_stoichiometry_float = stoic_r_float - stoic_f_float;
            stoichiometry_matrix[species_index, reaction_index] = net_stoichiometry

            ## Print N, Nf and Nr
            if net_stoichiometry!=0: 
                print "sbml.N(",species_index+1,",",reaction_index+1,") = ", net_stoichiometry, ";"
            if stoic_f!=0: 
                print "sbml.Nf(",species_index+1,",",reaction_index+1,") = ", stoic_f, ";"
            if stoic_r!=0: 
                print "sbml.Nr(",species_index+1,",",reaction_index+1,") = ", stoic_r, ";"

            if net_stoichiometry_float!=0: 
                print "sbml.N_float(",species_index+1,",",reaction_index+1,") = ", net_stoichiometry_float, ";"
            if stoic_f_float!=0: 
                print "sbml.Nf_float(",species_index+1,",",reaction_index+1,") = ", stoic_f_float, ";"
            if stoic_r_float!=0: 
                print "sbml.Nr_float(",species_index+1,",",reaction_index+1,") = ", stoic_r_float, ";"

    ## Species
    for species_ix, species_node in enumerate(sbml_model.getListOfSpecies()):
        #Name
        name = species_node.getName()
        print("sbml.species{{{0}}} = \"{1}\";".format(species_ix+1,name));
        
        # Units
        name = species_node.getSubstanceUnits()
        print("sbml.units{{{0}}} = \"{1}\";".format(species_ix+1,name));

        ## Initial conc
        name = species_node.getInitialConcentration()
        print("sbml.x_0{{{0}}} = \"{1}\";".format(species_ix+1,name));

    for reaction_ix, reaction_node in enumerate(sbml_model.getListOfReactions()):
        ## Name
        name = reaction_node.getName()
        print("sbml.reaction{{{0}}} = \"{1}\";".format(reaction_ix+1,name));

        ## KineticLaw
        law = reaction_node.getKineticLaw()
        formula = law.getFormula()
        print("sbml.kineticlaw{{{0}}} = \"{1}\";".format(reaction_ix+1,formula));
        # math = law.getMath()
        # print("sbml.kineticlaw{{{0}}} = \"{1}\";".format(reaction_ix+1,math));
        ##parameters = law.getListOfParameters()
        ##print("sbml.parameters{{{0}}} = \"{1}\";".format(reaction_ix+1,parameters));


    ## Indices of constant species
    print "sbml.I_const = ["
    for species_ix, species_node in enumerate(sbml_model.getListOfSpecies()):
        species_const = species_node.getConstant()
        if species_const==True:
            print("{0}".format(species_ix+1))
    print "];"

    return species, reactions, stoichiometry_matrix
                
if __name__ == '__main__':
    _main()
                                
                                

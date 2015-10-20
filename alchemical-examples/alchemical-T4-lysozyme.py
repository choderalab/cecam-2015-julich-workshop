#!/usr/bin/env python

"""
Create alchemical intermediates for default alchemical protocol for p-xylene in T4 lysozyme L99A in GBSA.

"""

from alchemy import AbsoluteAlchemicalFactory
from openmmtools import testsystems

# Create a reference system.
print "Creating a reference T4 lysozyme L99A system..."
complex = testsystems.LysozymeImplicit()
[reference_system, positions] = [complex.system, complex.positions]

# Create a factory to produce alchemical intermediates.
print "Creating an alchemical factory..."
receptor_atoms = range(0,2603) # T4 lysozyme L99A
ligand_atoms = range(2603,2621) # p-xylene
factory = AbsoluteAlchemicalFactory(reference_system, ligand_atoms=ligand_atoms)

# Get the default protocol for 'denihilating' in complex in explicit solvent.
protocol = factory.defaultComplexProtocolImplicit()

# Create the perturbed systems using this protocol.
print "Creating a perturbed system..."
systems = factory.createPerturbedSystems(protocol)
print "Done."

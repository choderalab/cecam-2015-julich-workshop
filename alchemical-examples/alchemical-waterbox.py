"""
Create alchemical intermediates for default alchemical protocol for one water in a water box.

"""
from alchemy import AbsoluteAlchemicalFactory, AlchemicalState
from openmmtools import testsystems

# Create a reference system.
print "Creating a water box..."
waterbox = testsystems.WaterBox()
[reference_system, positions] = [waterbox.system, waterbox.positions]

# Create a factory to produce alchemical intermediates.
print "Creating an alchemical factory..."
factory = AbsoluteAlchemicalFactory(reference_system, ligand_atoms=[0, 1, 2])

# Create a perturbed systems using this protocol.
print "Creating a perturbed system..."
alchemical_state = AlchemicalState()
alchemical_system = factory.createPerturbedSystem(alchemical_state)

# Perturb this system.
print "Perturbing the system..."
alchemical_state = AlchemicalState(lambda_sterics=0.90, lambda_electrostatics=0.90)
factory.perturbSystem(alchemical_system, alchemical_state)

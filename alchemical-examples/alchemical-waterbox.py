"""
Create alchemical intermediates for default alchemical protocol for one water in a water box.

"""
from alchemy import AbsoluteAlchemicalFactory
from openmmtools import testsystems

# Create a reference system.
print "Creating a water box..."
waterbox = testsystems.WaterBox()
[reference_system, positions] = [waterbox.system, waterbox.positions]

# Create a factory to produce alchemical intermediates.
print "Creating an alchemical factory..."
factory = AbsoluteAlchemicalFactory(reference_system, ligand_atoms=[0, 1, 2])

# Get the default protocol for 'denihilating' in solvent.
protocol = factory.defaultSolventProtocolExplicit()

# Create the perturbed systems using this protocol.
print "Creating alchemially perturbed systems..."
systems = factory.createPerturbedSystems(protocol)
print "Done."

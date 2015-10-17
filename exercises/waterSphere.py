from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

pdb = PDBFile('waterSphere.pdb')
forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
system = forcefield.createSystem(pdb.topology, nonbondedMethod=NoCutoff)
force = CustomExternalForce('10*max(0, r-1)^2; r=sqrt(x*x+y*y+z*z)')
for i in range(system.getNumParticles()):
    force.addParticle(i, ())
system.addForce(force)
integrator = LangevinIntegrator(1000*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)
simulation.reporters.append(PDBReporter('output.pdb', 100))
simulation.step(5000)

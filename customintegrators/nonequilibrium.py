#!/usr/bin/env python

"""
Nonequilibrium integrator based on Velocity Verlet implemented using CustomIntegrator

"""

from simtk import openmm, unit
from openmmtools import testsystems

timestep = 2.0 * unit.femtoseconds
temperature = 300.0 * unit.kelvin
nsteps = 100 # number of steps over which switching occurs
delta_K =  900.0 # change in spring constant K over nsteps (kJ/nm**2)

# Create a nonequilibrium integrator based on velocity Verlet.
integrator = openmm.CustomIntegrator(timestep)
integrator.addPerDofVariable("x1", 0) # holds positions prior to constraint correction
integrator.addComputeGlobal("testsystems_HarmonicOscillator_K", "testsystems_HarmonicOscillator_K + %f" % (delta_K/float(nsteps))) # update spring constant
integrator.addUpdateContextState() # allow other forces (e.g. MonteCarloBarostat) to update system
integrator.addComputePerDof("v", "v+0.5*dt*f/m") # velocity half kick
integrator.addComputePerDof("x", "x+dt*v") # position full kick
integrator.addComputePerDof("x1", "x") # store positions before position constraints are applied
integrator.addConstrainPositions() # apply position constraints (CCMA)
integrator.addComputePerDof("v", "v+0.5*dt*f/m") # velocity half kick
integrator.addComputePerDof("v", "v+(x-x1)/dt") # correct velocities for position constraints (RATTLE)
integrator.addConstrainVelocities() # constrain velocities (RATTLE)

# Run dynamics of a harmonic oscillator that changes spring constant dynamically
testsystem = testsystems.HarmonicOscillator()
context = openmm.Context(testsystem.system, integrator)
context.setPositions(testsystem.positions)
context.setVelocitiesToTemperature(temperature)
context.setParameter('testsystems_HarmonicOscillator_K', 100.0) # reset dynamic parameter that will be switched (kJ/nm**2)
state = context.getState(getEnergy=True)
initial_total_energy = state.getPotentialEnergy() + state.getKineticEnergy()
integrator.step(nsteps) # perform nonequilibrium switching
state = context.getState(getEnergy=True)
final_total_energy = state.getPotentialEnergy() + state.getKineticEnergy()
print "nonequilibrium work: %s" % str(final_total_energy - initial_total_energy)

#!/usr/bin/env python

"""
Velocity Verlet implemented using CustomIntegrator

"""

from simtk import openmm, unit
from openmmtools import testsystems

timestep = 1.0 * unit.femtoseconds

# Create a velocity Verlet integrator
integrator = openmm.CustomIntegrator(timestep)
integrator.addPerDofVariable("x1", 0) # holds positions prior to constraint correction
integrator.addUpdateContextState() # allow other forces (e.g. MonteCarloBarostat) to update system
integrator.addComputePerDof("v", "v+0.5*dt*f/m") # velocity half kick
integrator.addComputePerDof("x", "x+dt*v") # position full kick
integrator.addComputePerDof("x1", "x") # store positions before position constraints are applied
integrator.addConstrainPositions() # apply position constraints (CCMA)
integrator.addComputePerDof("v", "v+0.5*dt*f/m") # velocity half kick
integrator.addComputePerDof("v", "v+(x-x1)/dt") # correct velocities for position constraints (RATTLE)
integrator.addConstrainVelocities() # constrain velocities (RATTLE)

# Run dynamics
testsystem = testsystems.AlanineDipeptideImplicit()
context = openmm.Context(testsystem.system, integrator)
context.setPositions(testsystem.positions)
print "initial potential: " + str(context.getState(getEnergy=True).getPotentialEnergy())
integrator.step(10)
print "final potential:   " + str(context.getState(getEnergy=True).getPotentialEnergy())

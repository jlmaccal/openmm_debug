#!/usr/bin/env python

from simtk import openmm as mm
from simtk.openmm import app
from simtk import unit as u
import sys


RNG_SEED = 42


prm = app.AmberPrmtopFile('system.top')
crd = app.AmberInpcrdFile('system.mdcrd')

system = prm.createSystem(
    nonbondedMethod=app.CutoffPeriodic,
    nonbondedCutoff=1.0*u.nanometer,
    constraints=app.HBonds,
    rigidWater=True)

integrator = mm.LangevinIntegrator(300.*u.kelvin, 1.0/u.picosecond,
                                   1.0*u.femtosecond)
integrator.setRandomNumberSeed(RNG_SEED)

baro = mm.MonteCarloBarostat(1.0*u.atmosphere, 300.*u.kelvin, 25)
baro.setRandomNumberSeed(RNG_SEED)
system.addForce(baro)

platform = mm.Platform.getPlatformByName('CUDA')
properties = {'CudaPrecision': 'mixed'}

simulation = app.Simulation(prm.topology,
                            system,
                            integrator,
                            platform,
                            properties)

simulation.context.setPositions(crd.getPositions())
bv = crd.getBoxVectors()
simulation.context.setPeriodicBoxVectors(bv[0], bv[1], bv[2])

# print >>sys.stderr, 'minimizing'
# simulation.minimizeEnergy()
# print >>sys.stderr, 'done'

simulation.reporters.append(app.PDBReporter('traj.pdb', 1000))
simulation.reporters.append(app.StateDataReporter(
    sys.stderr, 1000, step=True, potentialEnergy=True,
    speed=True, separator='\t'))

for i in range(100):
    simulation.step(100)
    # state = simulation.context.getState(getPositions=True,
    #                                     getVelocities=True,
    #                                     getEnergy=True,
    #                                     enforcePeriodicBox=True)
    # simulation.context.setPositions(state.getPositions())
    # bv = state.getPeriodicBoxVectors()
    # simulation.context.setPeriodicBoxVectors(bv[0], bv[1], bv[2])

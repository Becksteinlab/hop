# $Id$
# specific simple testing code for density code

from hop.density import density_from_trajectory

def _setup_test():
    density = density_from_trajectory(psf="./NPT/ifabp_water.psf",
                                      dcd="./NPT/ifabp_water_1.dcd",
                                      delta=2.0,verbosity=3,
                                      metadata=dict(protein='I-FABP',species='water',
                                                    ligand=None,mutation=None,
                                                    PBC='cubic',temperature=300.0,
                                                    pressure=1.0,pdb='1IFC')
    )
    density.convert_density('TIP3P')

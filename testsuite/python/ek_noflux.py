import unittest as ut
import unittest_decorators as utx
import espressomd
import numpy as np
import espressomd.walberla
import espressomd.ekboundaries
import espressomd.shapes


@utx.skipIfMissingFeatures(["EK_WALBERLA"])
class EKNoFlux(ut.TestCase):
    BOX_L = 15.
    AGRID = 1.0
    DENSITY = 1
    DIFFUSION_COEFFICIENT = 0.1
    TIME = 50
    RADIUS = 5.

    system = espressomd.System(box_l=[BOX_L, BOX_L, BOX_L])
    system.time_step = 1.0
    system.cell_system.skin = 0.4

    def test_diffusion(self):
        """
        Testing the EK noflux boundaries to not leak density outside of a sphere.
        """

        espressomd.walberla.WalberlaBlockForest(
            box_size=self.system.box_l, ghost_layers=1, agrid=self.AGRID)

        ekspecies = espressomd.EKSpecies.EKSpecies(
            density=0.0, kT=0.0, diffusion=self.DIFFUSION_COEFFICIENT, valency=0.0)

        eksolver = espressomd.EKSpecies.EKNone()

        self.system.ekcontainer.add(ekspecies, tau=1.0, solver=eksolver)

        center = np.asarray(self.system.box_l / 2, dtype=np.int)

        ekspecies[center[0], center[1], center[2]].density = self.DENSITY

        sphere = espressomd.ekboundaries.EKBoundary(
            shape=espressomd.shapes.Sphere(
                center=self.system.box_l / 2,
                radius=self.RADIUS,
                direction=-1))
        self.system.ekboundaries.add(sphere)

        positions = np.empty((*self.system.box_l.astype(np.int), 3))
        positions[..., 2], positions[..., 1], positions[..., 0] = np.meshgrid(
            *map(lambda x: np.arange(0, x) - x / 2, self.system.box_l))
        positions += 0.5

        mask = np.linalg.norm(positions, axis=-1) < self.RADIUS

        self.system.integrator.run(self.TIME)

        simulated_density = np.copy(ekspecies[:, :, :].density)

        # check that the density is conserved
        np.testing.assert_almost_equal(
            np.sum(simulated_density), self.DENSITY, 10)
        if np.any(simulated_density < 0.):
            self.fail("EK density array contains negative densities!")

        # check that nothing leaked outside
        np.testing.assert_equal(simulated_density[~mask], 0.)


if __name__ == "__main__":
    ut.main()

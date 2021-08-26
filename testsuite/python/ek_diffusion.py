import unittest as ut
import unittest_decorators as utx
import espressomd
import numpy as np
import espressomd.walberla
import scipy.optimize


@utx.skipIfMissingFeatures(["EK_WALBERLA"])
class EKDiffusion(ut.TestCase):
    BOX_L = 31.
    AGRID = 1.0
    DENSITY = 1
    DIFFUSION_COEFFICIENT = 0.1
    TIME = 150

    system = espressomd.System(box_l=[BOX_L, BOX_L, BOX_L])
    system.time_step = 1.0
    system.cell_system.skin = 0.4

    def analytical_density(self, pos: np.ndarray, time: int, D: float):
        return (4 * np.pi * D * time)**(-3 / 2) * \
            np.exp(-np.sum(np.square(pos), axis=-1) / (4 * D * time))

    def test_diffusion(self):
        """
        Testing EK for simple diffusion of a point droplet
        """

        espressomd.walberla.WalberlaBlockForest(
            box_size=self.system.box_l, ghost_layers=1, agrid=self.AGRID)

        ekspecies = espressomd.EKSpecies.EKSpecies(
            density=0.0, kT=0.0, diffusion=self.DIFFUSION_COEFFICIENT, valency=0.0, advection=False, friction_coupling=False)

        eksolver = espressomd.EKSpecies.EKNone()

        self.system.ekcontainer.add(ekspecies, tau=1.0, solver=eksolver)

        center = np.asarray(self.system.box_l / 2, dtype=np.int)

        ekspecies[center].density = self.DENSITY

        # check that the density in the domain is what is expected
        np.testing.assert_almost_equal(
            np.sum(ekspecies[:, :, :].density), self.DENSITY, 10)

        # TODO: replace that when the blockforest is able to return the
        # dimensions
        positions = np.empty((*self.system.box_l.astype(np.int), 3))
        positions[..., 2], positions[..., 1], positions[..., 0] = np.meshgrid(
            *map(lambda x: np.arange(0, x) - x / 2, self.system.box_l))
        positions += 0.5

        self.system.integrator.run(self.TIME)

        simulated_density = np.copy(ekspecies[:, :, :].density)

        # check that the density is conserved
        np.testing.assert_almost_equal(
            np.sum(simulated_density), self.DENSITY, 10)
        if np.any(simulated_density < 0.):
            self.fail("EK density array contains negative densities!")

        # check that the maximum is in the right place
        peak = np.unravel_index(
            np.argmax(
                simulated_density,
                axis=None),
            self.system.box_l.astype(
                np.int))
        np.testing.assert_equal(peak, self.system.box_l / 2 - 0.5)

        calc_density = self.analytical_density(
            positions, self.TIME, self.DIFFUSION_COEFFICIENT)
        target = [self.TIME, self.DIFFUSION_COEFFICIENT]

        popt, _ = scipy.optimize.curve_fit(self.analytical_density,
                                           positions.reshape(-1, 3),
                                           simulated_density.reshape(-1),
                                           p0=target,
                                           bounds=([0, 0], [np.inf, np.inf]))

        np.testing.assert_allclose(popt[0], self.TIME, rtol=0.1)
        np.testing.assert_allclose(
            popt[1], self.DIFFUSION_COEFFICIENT, rtol=0.1)
        np.testing.assert_allclose(
            calc_density,
            simulated_density,
            atol=1e-5,
            rtol=0)


if __name__ == "__main__":
    ut.main()

# Measuring the force on a single sphere immersed in a fluid with
# fixed velocity boundary conditions created by two
# walls at finite distance.
# The force is compared to th analytical result F=6 pi eta r v
# i.e. the stokes force on the particles.


# We create a box of size box_width x box_width x box_length and
# place an object in the center. We measure the drag force
# in z direction. We create walls in the xz and yz plane at the box
# boundaries, where the velocity is fixed to $v.
#
from espressomd import System, lb, lbboundaries, shapes, has_features
import unittest as ut
import numpy as np
import sys


@ut.skipIf(not has_features(["LB", "LB_BOUNDARIES"]),
           "Features not available, skipping test!")
class Stokes(ut.TestCase):

    es = System()

    def run_stokes(self,lb_class):
        # System setup
        #system = System()
        system = self.es
        agrid = 1
        radius = 5.5
        box_width = 44 
        real_width = box_width + 2 * agrid
        box_length = 44 
        system.box_l = [real_width, real_width, box_length]
        system.time_step = 0.5
        system.cell_system.skin = 0.4

        # The temperature is zero.
        system.thermostat.set_lb(kT=0)

        # LB Parameters
        v = [0, 0, 0.01]  # The boundary slip
        kinematic_visc = 2.5

        # Invoke LB fluid
        lbf = lb_class(visc=kinematic_visc, dens=1,
                             agrid=agrid, tau=system.time_step, fric=1)
        system.actors.add(lbf)

        # Setup walls
        walls = [None] * 4
        walls[0] = lbboundaries.LBBoundary(shape=shapes.Wall(
            normal=[-1, 0, 0], dist=-(1 + box_width)), velocity=v)
        walls[1] = lbboundaries.LBBoundary(
            shape=shapes.Wall(
                normal=[
                    1,
                    0,
                    0],
                dist=1),
            velocity=v)
        walls[2] = lbboundaries.LBBoundary(shape=shapes.Wall(
            normal=[0, -1, 0], dist=-(1 + box_width)), velocity=v)
        walls[3] = lbboundaries.LBBoundary(
            shape=shapes.Wall(
                normal=[
                    0,
                    1,
                    0],
                dist=1),
            velocity=v)

        for wall in walls:
            system.lbboundaries.add(wall)

        # setup sphere without slip in the middle
        sphere = lbboundaries.LBBoundary(shape=shapes.Sphere(
            radius=radius, center=[real_width / 2] * 2 + [box_length / 2], direction=1))

        system.lbboundaries.add(sphere)

        def size(vector):
            tmp = 0
            for k in vector:
                tmp += k * k
            return np.sqrt(tmp)

        stokes_force = 6 * np.pi * kinematic_visc * radius * size(v)
        print("Stokes' Law says: f=%f" % stokes_force)
        system.integrator.run(800)
        first_force=size(sphere.get_force())
        system.integrator.run(120)
        force = sphere.get_force()
        print("Measured force: f=%f" % size(force))
        second_force=size(force)
        
        # Tidy before asserts
        for wall in walls:
            system.lbboundaries.add(wall)
        # Lb fluid deactivation code missing. Un-comment once implemented
        #del system.actors[0]
        system.part.clear()

        self.assertAlmostEqual(first_force,second_force,places=5)


        print("\nIntegration finished.")

        # get force that is exerted on the sphere
        self.assertLess(abs(1.0 - size(force) / stokes_force), 0.06)

    def test_cpu(self):
        self.run_stokes(lb.LBFluid)
    


if __name__ == "__main__":
    ut.main()

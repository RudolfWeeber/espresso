#
# Copyright (C) 2026 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
"""
Regression test for an MPI deadlock in Stokesian Dynamics.

When a particle has a type for which no radius was registered, the core
used to call ``radii.at(p.type)`` only inside the rank-0-only block of
``StokesianDynamics::propagate_vel_pos``, *after* the collective
``gather_buffer`` but *before* the collective ``scatter_buffer``. The
resulting ``std::out_of_range`` unwound out of the integration loop on
rank 0, while the other ranks blocked forever in the matching
``MPI_Scatterv`` -> indefinite deadlock.

This test is registered ``NO_MPI`` so that the parent process is a plain
single process. It then launches the offending scenario as a child
``mpiexec -n 2`` job under a hard timeout. OpenMPI forbids recursive
``mpiexec`` calls, hence the deadlock cannot be reproduced from within an
already-MPI parent. A hang (subprocess timeout) is the failure signature
of the bug; a clean coordinated runtime error (the child exits and prints
``OK``) is the fixed behaviour.
"""
import unittest as ut
import unittest_decorators as utx
import os
import sys
import shutil
import pathlib
import subprocess


@utx.skipIfMissingFeatures("STOKESIAN_DYNAMICS")
class StokesianMissingRadius(ut.TestCase):

    # pypresso exports PYTHONPATH, so the plain interpreter in `sys.executable`
    # can import espressomd in the child (same pattern as `caliper.py`)
    interpreter = sys.executable
    child = str(pathlib.Path(__file__).resolve().parent /
                "stokesian_missing_radius_child.py")

    @staticmethod
    def _clean_mpi_env():
        # importing espressomd already initialized MPI in this (parent)
        # process; OpenMPI then refuses the nested ``mpiexec`` ("does not
        # support recursive calls") unless the inherited MPI environment is
        # stripped from the child's environment.
        return {k: v for k, v in os.environ.items()
                if not (k.startswith("OMPI_") or k.startswith("PMIX_") or
                        k.startswith("PMI_"))}

    def test_missing_radius_no_deadlock(self):
        """
        A particle whose type has no registered radius must trigger a
        coordinated error on all ranks, not an MPI deadlock.
        """
        mpiexec = shutil.which("mpiexec")
        if mpiexec is None:
            self.skipTest("mpiexec not available")

        cmd = [mpiexec, "-n", "2", "--oversubscribe", "--bind-to", "none",
               self.interpreter, self.child]
        try:
            result = subprocess.run(
                cmd, timeout=45, capture_output=True, text=True,
                env=self._clean_mpi_env())
        except subprocess.TimeoutExpired as err:
            self.fail(
                "Stokesian Dynamics deadlocked on a particle with a missing "
                f"radius (subprocess timed out).\nstdout:\n{err.stdout}\n"
                f"stderr:\n{err.stderr}")

        combined = (result.stdout or "") + (result.stderr or "")
        # The child must NOT have completed integration silently and must NOT
        # have aborted with a raw, uncoordinated std::out_of_range.
        self.assertNotIn("unordered_map::at", combined,
                         msg=f"raw uncoordinated abort:\n{combined}")
        self.assertNotIn("NO ERROR RAISED", combined,
                         msg=f"missing radius went unnoticed:\n{combined}")
        # The child reports OK once it sees a clean, coordinated runtime error
        # mentioning the missing radius.
        self.assertIn("OK: coordinated runtime error raised", combined,
                      msg=f"unexpected child output:\n{combined}")


if __name__ == "__main__":
    ut.main()

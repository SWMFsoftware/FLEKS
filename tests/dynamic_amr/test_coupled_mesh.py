import pathlib
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tests.dynamic_amr.check_coupled_mesh import check_mesh_change


class CoupledMeshTest(unittest.TestCase):
    def test_second_session_changes_fine_grid(self):
        log = """----- Starting Session            1  ------
FLEKS0: Domain::regrid is called
FLEKS0 fi:  iLev = 1  # of boxes = 144  # of cells = 9216
----- Starting Session            2  ------
FLEKS0: Domain::regrid is called
FLEKS0 fi:  iLev = 1  # of boxes = 64  # of cells = 4096
"""
        self.assertEqual(check_mesh_change(log), (144, 9216, 64, 4096))

    def test_same_fine_grid_fails(self):
        log = """----- Starting Session 1 ------
FLEKS0: Domain::regrid is called
FLEKS0 fi: iLev = 1 # of boxes = 64 # of cells = 4096
----- Starting Session 2 ------
FLEKS0: Domain::regrid is called
FLEKS0 fi: iLev = 1 # of boxes = 64 # of cells = 4096
"""
        with self.assertRaisesRegex(ValueError, "did not change"):
            check_mesh_change(log)

    def test_missing_second_regrid_fails(self):
        log = """----- Starting Session 1 ------
FLEKS0: Domain::regrid is called
FLEKS0 fi: iLev = 1 # of boxes = 64 # of cells = 4096
----- Starting Session 2 ------
FLEKS0 fi: iLev = 1 # of boxes = 32 # of cells = 2048
"""
        with self.assertRaisesRegex(ValueError, "regrid"):
            check_mesh_change(log)


if __name__ == "__main__":
    unittest.main()

"""Regression for immutable named regions across parameter sessions."""

import pathlib
import subprocess
import tempfile
import unittest


ROOT = pathlib.Path(__file__).resolve().parents[2]
EXECUTABLE = ROOT / "bin" / "FLEKS.exe"
PARAM = pathlib.Path(__file__).with_name("PARAM.in")


class RegionNameTest(unittest.TestCase):
    def test_new_shape_waits_for_refine_region_selection(self):
        first_session = PARAM.read_text().split("#RUN", 1)[0]
        first_session = first_session.replace(
            "5                       MaxIter", "0                       MaxIter"
        )
        later_sessions = """#REGION
right                   name
box                     shape
0.0                     min
1000.0                  max
-10.0                   min
10.0                    max
-10.0                   min
10.0                    max

#STOP
1                       MaxIter
-1.0                    TimeMax

#RUN

#REFINEREGION
0                       iLev
+right                  regions

#STOP
2                       MaxIter
-1.0                    TimeMax

#RUN

#REFINEREGION
0                       iLev
none                    regions

#STOP
3                       MaxIter
-1.0                    TimeMax

#RUN
"""
        with tempfile.TemporaryDirectory() as workdir:
            run_dir = pathlib.Path(workdir)
            (run_dir / "PARAM.in").write_text(
                first_session + "#RUN\n\n" + later_sessions
            )
            result = subprocess.run(
                [str(EXECUTABLE)], cwd=run_dir,
                capture_output=True, text=True, timeout=120,
            )
        output = result.stdout + result.stderr
        self.assertEqual(result.returncode, 0, output[-4000:])
        # Initial grid, switch to right, and clear refinement. Defining right
        # by itself must not cause an extra regrid.
        self.assertEqual(output.count("Domain::regrid is called"), 3)

    def test_redefinition_in_later_session_is_rejected(self):
        first_session = PARAM.read_text().split("#RUN", 1)[0]
        first_session = first_session.replace(
            "5                       MaxIter", "0                       MaxIter"
        )

        for x_min, x_max in [(-1000.0, 0.0), (0.0, 1000.0)]:
            with self.subTest(x_min=x_min, x_max=x_max):
                second_session = f"""#REGION
left                    name
box                     shape
{x_min}                  min
{x_max}                  max
-10.0                   min
 10.0                   max
-10.0                   min
 10.0                   max

#STOP
0                       MaxIter
-1.0                    TimeMax

#RUN
"""
                with tempfile.TemporaryDirectory() as workdir:
                    run_dir = pathlib.Path(workdir)
                    (run_dir / "PARAM.in").write_text(
                        first_session + "#RUN\n\n" + second_session
                    )
                    result = subprocess.run(
                        [str(EXECUTABLE)], cwd=run_dir,
                        capture_output=True, text=True, timeout=120,
                    )
                output = result.stdout + result.stderr
                self.assertNotEqual(result.returncode, 0, output)
                self.assertIn("Duplicate #REGION name", output)


if __name__ == "__main__":
    unittest.main()

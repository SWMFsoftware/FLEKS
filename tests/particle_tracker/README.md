# Particle Tracker Tests

Standalone test suite for pure particle tracking (`PIC F` with `SOLVEEM F`) and
test-particle pusher/output infrastructure.

## Variants

1. `PARAM.in` — Relativistic Boris (`TPPUSHER boris`) crossed `E x B` drift benchmark
   with two proton species (force-free `0.8 c` beam and cycloid from rest). Verifies
   log suppression (`log_pic_n*.log` omitted), standard 13-column `PBE` binary output,
   and `flekspy.FLEKSTP` compatibility.
2. `PARAM.in.vay` — Relativistic Vay (2008) pusher (`TPPUSHER vay`) on the same setup,
   verifying exact force-free `E + v x B = 0` cancellation and boosted-frame energy
   conservation.
3. `PARAM.in.hc` — Relativistic Higuera-Cary (2017) pusher (`TPPUSHER hc`) with runtime
   `TPOUTPUT` (`PBE vexb energy gamma`, 18 columns) flushed across multiple windows
   (`nRecordFlush = 8`).
4. `PARAM.in.amr` — Two-level 3D AMR test with static `DIPOLE`, absorbing `BODY`,
   specularly reflecting X walls (`PARTICLEBOXBOUNDARY reflect`), `E = 0` everywhere, and
   all derived diagnostics (`PBE drifts vexb mu pitch energy gamma`, 26 columns) verified
   by `validate.py`.

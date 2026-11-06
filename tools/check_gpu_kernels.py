#!/usr/bin/env python3
"""Static checks that keep FLEKS kernels GPU-safe.

FLEKS runs on CPUs and on NVIDIA GPUs from one source tree. A kernel passed to
amrex::ParallelFor (and friends) runs on the device in a GPU build, so it must

  1. be a device lambda: ``[=] AMREX_GPU_DEVICE (...)``. A host lambda such as
     ``[&](int i, int j, int k)`` does not compile for the GPU, and a
     by-reference capture would dereference host stack memory on the device;
  2. not capture ``this``. Reading a class member inside a ``[=]`` lambda
     captures ``this``, which points to host memory. nvcc rejects it
     (--Werror ext-lambda-captures-this) but gcc does not, so this script
     cannot see it; the CUDA build in CI catches it.

GPU-specific code (``#ifdef AMREX_USE_GPU`` / ``AMREX_USE_CUDA``) is allowed
only in the files listed in GPU_BACKEND_FILES so that the CPU and the GPU
builds keep one code path.

Usage: python3 tools/check_gpu_kernels.py [files...]   (default: include, src)
Exit status is 1 if any problem is found.
"""

import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Launch functions whose lambda argument runs on the device in a GPU build.
LAUNCHERS = (
    "ParallelFor",
    "ParallelForRNG",
    "ParallelForOMP",
    "ParticleReduce",
    "ParticleToMesh",
    "MeshToParticle",
    "ParReduce",
    "AMREX_HOST_DEVICE_PARALLEL_FOR_1D",
    "AMREX_HOST_DEVICE_PARALLEL_FOR_3D",
    "AMREX_HOST_DEVICE_PARALLEL_FOR_4D",
)

# Files that may contain GPU-only preprocessor branches.
GPU_BACKEND_FILES = {
    "include/GpuBackend.h",
}

LAUNCH_RE = re.compile(r"\b(?:amrex::)?(" + "|".join(LAUNCHERS) +
                       r")\s*\(|\.(eval)\s*\(")  # .eval: amrex::ReduceOps
LAMBDA_RE = re.compile(
    r"\[([^\[\]]*)\]\s*(AMREX_GPU_DEVICE|AMREX_GPU_HOST_DEVICE)?\s*\(")
GPU_IFDEF_RE = re.compile(
    r"^\s*#\s*(if|ifdef|ifndef|elif)\b.*\bAMREX_USE_(GPU|CUDA)\b")


def strip_comments(text):
    """Blank out comments and string literals, keeping line numbers."""
    out = []
    i, n = 0, len(text)
    while i < n:
        c = text[i]
        if text.startswith("//", i):
            j = text.find("\n", i)
            j = n if j < 0 else j
            out.append(" " * (j - i))
            i = j
        elif text.startswith("/*", i):
            j = text.find("*/", i + 2)
            j = n if j < 0 else j + 2
            out.append(re.sub(r"[^\n]", " ", text[i:j]))
            i = j
        elif c in "\"'":
            j = i + 1
            while j < n and text[j] != c:
                j += 2 if text[j] == "\\" else 1
            j = min(j + 1, n)
            out.append(c + " " * (j - i - 2) + c if j - i >= 2 else text[i:j])
            i = j
        else:
            out.append(c)
            i += 1
    return "".join(out)


def call_arguments(text, start):
    """Return the text between the parenthesis at ``start`` and its match."""
    depth = 0
    for j in range(start, len(text)):
        if text[j] == "(":
            depth += 1
        elif text[j] == ")":
            depth -= 1
            if depth == 0:
                return text[start + 1:j]
    return text[start + 1:]


def check_file(path):
    rel = os.path.relpath(path, ROOT)
    with open(path, encoding="utf-8", errors="replace") as f:
        raw = f.read()
    text = strip_comments(raw)
    problems = []

    for m in LAUNCH_RE.finditer(text):
        # interp.ParticleToMesh(...) is a method of a particle interpolator
        # called inside a kernel, not a launch.
        before = text[max(0, m.start() - 2):m.start()]
        if m.group(1) and (before.endswith(".") or before.endswith("->")):
            continue
        args = call_arguments(text, m.end() - 1)
        lam = LAMBDA_RE.search(args)
        if lam is None:
            continue  # functor or named lambda; nvcc checks those
        line = text.count("\n", 0, m.start()) + 1
        name = m.group(1) or m.group(2)
        capture, qualifier = lam.group(1).strip(), lam.group(2)
        if qualifier is None:
            problems.append(f"{rel}:{line}: {name} lambda is not "
                            "AMREX_GPU_DEVICE")
        elif "&" in capture:
            problems.append(f"{rel}:{line}: {name} device lambda "
                            f"captures by reference [{capture}]")

    if rel not in GPU_BACKEND_FILES:
        for line_no, line in enumerate(text.splitlines(), 1):
            if GPU_IFDEF_RE.match(line):
                problems.append(f"{rel}:{line_no}: GPU-only preprocessor "
                                "branch outside the GPU backend files")
    return problems


def main(argv):
    if len(argv) > 1:
        files = argv[1:]
    else:
        files = []
        for sub in ("include", "src"):
            for dirpath, _, names in os.walk(os.path.join(ROOT, sub)):
                files += [os.path.join(dirpath, f) for f in sorted(names)
                          if f.endswith((".h", ".cpp", ".H"))]
    problems = []
    for path in files:
        problems += check_file(path)
    for p in problems:
        print(p)
    if problems:
        print(f"\n{len(problems)} GPU kernel problem(s) found. See the "
              "docstring of tools/check_gpu_kernels.py.")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

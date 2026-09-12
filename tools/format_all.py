import os
import re
import shutil
import subprocess
import sys

EXPECTED_CLANG_FORMAT_VERSION = "20.1.8"

def check_formatters():
    clang_format_bin = shutil.which("clang-format")
    if not clang_format_bin:
        print(
            f"WARNING: 'clang-format' not found in PATH.\n"
            f"         C++ files will not be formatted!\n"
            f"         Install CI-matching version via:\n"
            f"           python3 -m pip install clang-format=={EXPECTED_CLANG_FORMAT_VERSION}\n"
            f"         (or: uv tool install clang-format=={EXPECTED_CLANG_FORMAT_VERSION})\n",
            file=sys.stderr,
        )
    else:
        try:
            res = subprocess.run(
                [clang_format_bin, "--version"],
                capture_output=True,
                text=True,
                check=True,
            )
            match = re.search(r"version\s+([0-9]+\.[0-9]+\.[0-9]+)", res.stdout)
            if match:
                version = match.group(1)
                if version != EXPECTED_CLANG_FORMAT_VERSION:
                    print(
                        f"WARNING: clang-format version is {version}, but CI uses {EXPECTED_CLANG_FORMAT_VERSION}.\n"
                        f"         Formatting differences may occur in CI.\n"
                        f"         Install matching version via: python3 -m pip install clang-format=={EXPECTED_CLANG_FORMAT_VERSION}\n",
                        file=sys.stderr,
                    )
        except Exception:
            pass

    findent_bin = shutil.which("findent")
    if not findent_bin:
        print(
            "NOTICE: 'findent' not found in PATH. Fortran files will not be formatted.\n"
            "        Install via: python3 -m pip install findent\n"
            "        (or: uv tool install findent)\n",
            file=sys.stderr,
        )

def process_file(file_path):
    print(f"Processing {file_path}")
    if file_path.endswith((".f90", ".F90")):
        # run findent
        temp_file = file_path + ".tmp"
        with open(file_path, "r") as f_in, open(temp_file, "w") as f_out:
            args = ["findent", "-i3", "-r2", "-m2", "-k5", "-c3", "-C2", "-j2", "-a2"]
            try:
                subprocess.run(args, stdin=f_in, stdout=f_out, check=True)
            except FileNotFoundError:
                pass
        if os.path.exists(temp_file) and os.path.getsize(temp_file) > 0:
            os.replace(temp_file, file_path)
        else:
            if os.path.exists(temp_file):
                os.remove(temp_file)

    elif file_path.endswith((".cpp", ".h")):
        try:
            subprocess.run(["clang-format", "-i", file_path], check=True)
        except FileNotFoundError:
            pass

    # for all files, remove trailing whitespace and ensure newline at EOF
    with open(file_path, "r", encoding="utf-8") as f:
        lines = f.readlines()

    # Strip trailing whitespaces
    lines = [line.rstrip() + "\n" for line in lines]

    # Ensure file ends with exactly one newline if not empty
    if lines:
        while len(lines) > 1 and lines[-1] == "\n" and lines[-2] == "\n":
            lines.pop()

    with open(file_path, "w", encoding="utf-8") as f:
        f.writelines(lines)

def collect_files(targets):
    valid_exts = (".h", ".cpp", ".f90", ".F90")
    files_to_format = []
    for t in targets:
        if os.path.isfile(t):
            if t.endswith(valid_exts):
                files_to_format.append(t)
        elif os.path.isdir(t):
            for root, _, files in os.walk(t):
                for file in files:
                    if file.endswith(valid_exts):
                        files_to_format.append(os.path.join(root, file))
        else:
            print(f"Target not found: {t}", file=sys.stderr)
    return files_to_format

def main():
    check_formatters()

    targets = sys.argv[1:]
    if targets:
        files = collect_files(targets)
    else:
        repo_dirs = ["src", "include", "srcInterface"]
        files = collect_files(repo_dirs)

    for file_path in files:
        process_file(file_path)

if __name__ == "__main__":
    main()


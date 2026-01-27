import json, os, glob, subprocess, shutil

def generate():
    root = os.getcwd()
    swmf = os.path.abspath(os.path.join(root, "../../"))
    compiler = shutil.which("mpicxx") or "/opt/local/bin/mpicxx"

    # Get MPI flags
    mpi_out = subprocess.check_output([compiler, "-show"], text=True)
    includes = [p for p in mpi_out.split() if p.startswith("-I")]
    
    # Get COMPONENT flag
    comp_out = subprocess.check_output(["grep", "^COMPONENT", os.path.join(root, "Makefile.def")], text=True).strip()
    flags = [f"-D_{comp_out.split('=')[1].strip()}_COMPONENT_"]

    includes += ["-I../include", 
                 f"-I{os.path.join(swmf, 'share/Library/src')}",
                 f"-I{os.path.join(swmf, 'util/AMREX/InstallDir/include')}"]

    base_cmd = " ".join([compiler] + includes + flags)
    src_dir = os.path.join(root, "src")
    
    entries = [{
        "directory": src_dir,
        "command": f"{base_cmd} -c {os.path.basename(f)} -o {os.path.splitext(os.path.basename(f))[0]}.o",
        "file": os.path.basename(f)
    } for f in sorted(glob.glob(os.path.join(src_dir, "*.cpp")))]

    with open(os.path.join(root, "compile_commands.json"), 'w') as f:
        json.dump(entries, f, indent=4)
        
    print(f"Generated compile_commands.json with {len(entries)} entries.")

if __name__ == "__main__":
    generate()

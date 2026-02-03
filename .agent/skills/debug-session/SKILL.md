---
name: Debug Session
description: Set up and launch a debugging session for FLEKS using VS Code or lldb
---

# Debug Session

This skill helps set up and run debugging sessions for FLEKS.

## Prerequisites

- FLEKS built with debug symbols (`-g` flag)
- For VS Code: CodeLLDB or C/C++ extension installed
- For command line: `lldb` (macOS) or `gdb` (Linux)

## VS Code Debugging

### 1. Ensure Debug Build

The project should be compiled with debug flags. Check `Makefile.conf` for optimization level.

### 2. Launch Configuration

The existing launch configuration is in `.vscode/launch.json`:
- **Program**: `${workspaceFolder}/../../run/SWMF.exe`
- **Working Directory**: `${workspaceFolder}/../../run/`
- **Debugger**: `lldb` on macOS, `gdb` on Linux

### 3. Start Debugging

1. Open VS Code
2. Press `F5` or go to Run → Start Debugging
3. Select "cpp build and debug active file"

### 4. Setting Breakpoints

- Click in the gutter (left margin) next to line numbers
- Or use `Debug: Toggle Breakpoint` command
- Conditional breakpoints: Right-click breakpoint → Edit Breakpoint

## Command Line Debugging (lldb)

### 1. Start lldb

```bash
cd /Users/yuxichen/shock/SWMF/run
lldb ./SWMF.exe
```

### 2. Common lldb Commands

| Command | Description |
|---------|-------------|
| `b <file>:<line>` | Set breakpoint at file:line |
| `b <function>` | Set breakpoint at function |
| `r` | Run program |
| `c` | Continue execution |
| `n` | Step over (next) |
| `s` | Step into |
| `finish` | Step out of function |
| `p <var>` | Print variable value |
| `bt` | Show backtrace |
| `frame select <n>` | Select stack frame |
| `quit` | Exit debugger |

### 3. Example Session

```bash
lldb ./SWMF.exe
(lldb) b Pic.cpp:100
(lldb) r
(lldb) p variableName
(lldb) bt
(lldb) c
```

## MPI Debugging

For debugging parallel runs, you have several options:

### Option 1: Attach to Running Process

```bash
# Run MPI program
mpirun -np 4 ./SWMF.exe &

# Find process ID
ps aux | grep SWMF

# Attach lldb
lldb -p <pid>
```

### Option 2: Debug Single Rank

```bash
mpirun -np 4 xterm -e lldb ./SWMF.exe
```

## Useful Watch Expressions

For physics simulations, useful things to watch:
- Time step counters
- Grid dimensions
- Particle counts
- Field values at specific indices

## Troubleshooting

| Issue | Solution |
|-------|----------|
| No symbols | Rebuild with `-g` flag |
| Can't find source | Check source path in launch.json |
| Breakpoint not hit | Verify code path is executed |
| Optimized away | Reduce optimization level |

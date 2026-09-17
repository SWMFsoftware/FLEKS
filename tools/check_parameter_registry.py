#!/usr/bin/env python3
"""Verify registered PARAM.in commands are documented in PARAM.XML."""

import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def registered_commands() -> set[str]:
    domain = (ROOT / "src" / "Domain.cpp").read_text(encoding="utf-8")
    source = (ROOT / "userfiles" / "ExoSource.h").read_text(encoding="utf-8")
    commands = set(re.findall(r'\{"(#[A-Z0-9]+)",\s*ParameterOwner::', domain))
    commands.update(re.findall(r'commands\.push_back\("(#[A-Z0-9]+)"\)', source))
    return commands


def documented_commands() -> set[str]:
    text = (ROOT / "PARAM.XML").read_text(encoding="utf-8")
    commands = {
        command.upper()
        for command in re.findall(r"#[A-Z][A-Z0-9_]*", text.upper())
    }
    for tag in re.findall(r"<command\b[^>]*>", text, flags=re.IGNORECASE):
        for value in re.findall(r'(?:name|alias)="([^"]*)"', tag):
            for name in value.split(","):
                name = name.strip().upper()
                if name:
                    commands.add("#" + name.lstrip("#"))
    return commands


def main() -> int:
    registered = registered_commands()
    documented = documented_commands()
    missing = sorted(registered - documented)
    if missing:
        print("Registered commands missing from PARAM.XML:")
        print("\n".join(f"  {command}" for command in missing))
        return 1
    print(f"Verified {len(registered)} registered commands against PARAM.XML.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

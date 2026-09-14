# FLEKS Parameter Reference

> Canonical: `PARAM.XML`, rendered by `make PDF` into `doc/USERMANUAL.pdf`.
> This page is only a navigation aid — read the `<command>` block in `PARAM.XML`
> before using or changing a command.

## Major command groups

| Group | Commands |
|---|---|
| Output | `#SAVEPLOT`, `#MONITOR`, `#SAVELOG`, `#NOUTFILE` |
| Scheme | `#PIC`, `#TIMESTEPPING`, `#DISCRETIZATION`, `#EFIELDSOLVER`, `#DIVE`, `#DIVB` |
| Hybrid PIC | `#HYBRIDPIC`, `#RESISTIVITY`, `#HALLTERM`, `#ELECTRONTEMPERATURE`, `#HYPERRESISTIVITY`, `#BSUBCYCLE`, `#MINIMUMDENSITY`, `#FIELDINTEGRATOR`, `#AVGFIELDB`, `#SMOOTHMOMENTS` |
| Particles | `#PARTICLES`, `#RESAMPLING`, `#FASTMERGE`, `#VACUUM`, `#PARTICLETRACKER` |
| Initial / boundary | `#GEOMETRY`, `#NCELL`, `#REGION`, `#BC` |
| Coupling | `#OHION`, `#CHARGEEXCHANGE`, `#MAXCHARGEEXCHANGERATE` |

## Conventions

- Per-domain overrides use `multiple="T"` plus `COMMANDNAME_FLEKS{0,1,2}`
  aliases.
- Do not use `#` to reference a command inside a `PARAM.in` comment — it is
  parsed as a command.
- Adding a new command end-to-end: `.agents/workflows/add-param.md`.

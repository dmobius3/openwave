# M8.8 Consulted Files

Files read during this cleanroom session, per protocol § 8 step 9.

## Input packets (within cleanroom)

1. **m8_5a_packet.json** — Group packet defining the binary icosahedral group 2I.
   Contains: coefficient field Q(φ), two quaternion generators, quaternion basis.
   Used for: element enumeration, SU(2) embedding, representation construction.

2. **m8_8_construction_packet.json** — Construction packet with based chain complex.
   Contains: abstract generators (s=118, t=80), free ranks [1,2,2,1],
   boundary maps d1 (2×1), d2 (2×2), d3 (1×2) over Z[2I], basing conventions.
   Used for: chain complex validation, boundary map evaluation, torsion computation.

## Protocol and task files (within cleanroom)

3. **TASK.md** — Operational instructions for the cleanroom reproduction.

4. **PROTOCOL.md** — M8.8 Independent-Method Reproduction Protocol.
   Sections consulted: § 4.2 (enumeration), § 5.4 (orientation involution),
   § 5.5 (output schema), § 8 (adjudication sequence).

## Files NOT consulted

- No files outside `/Users/blake/Desktop/OpenWave/M8.8_CLEANROOM/`
- No M8.3 implementation files
- No spectral data, heat-kernel data, or zeta function values
- No published torsion tables or ratio tables for 2I
- No mode-identity-theory artifact
- No answer packet

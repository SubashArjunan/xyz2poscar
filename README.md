# xyz2poscar

A small, dependency-free Bash/awk script that converts an `.xyz` file into a VASP `POSCAR`
with **Selective Dynamics**, so you can freeze chosen atoms (or whole elements) in one command.
It also takes care of the usual chores when preparing a molecule or cluster for VASP
(including MD runs): moving the structure into the positive octant, building a box, and
centring the molecule with a chosen amount of vacuum.

## Features

- Freeze atoms by index or range (`-f "1-5,10"`), by element (`-e "Au,O"`), or both (union)
- Frozen atoms get `F F F`, all others `T T T`
- Atoms are grouped by element as VASP requires; the original `.xyz` index is kept as a
  comment on every line (`! atom 7 (Au)`), so `-f` indices always refer to the input file
- Negative coordinates are rigidly shifted into the +x, +y, +z octant (geometry unchanged)
- Cell handling, in order of priority:
  1. an explicit cell with `-c`
  2. an extended-XYZ `Lattice="..."` header (as written by ASE and many MD tools)
  3. an automatic orthorhombic box
- `-v VAC` builds a box with `VAC` angstrom of vacuum on each side and centres the molecule
- `-m` centres the molecule in any cell
- Cartesian (default) or Direct (fractional) output

## Requirements

- Bash 3.2 or newer
- A POSIX `awk` (tested with mawk 1.3.4 and GNU awk)

No Python or other dependencies.

## Installation

```bash
git clone https://github.com/<your-user>/xyz2poscar.git
cd xyz2poscar
chmod +x xyz2poscar.sh        # only if the executable bit was lost
# optional: put it on your PATH
cp xyz2poscar.sh ~/bin/xyz2poscar
```

## Usage

```bash
./xyz2poscar.sh -i structure.xyz [-c "ax ay az bx by bz cx cy cz"] \
                [-o POSCAR] [-t Cartesian|Direct] \
                [-f "1-5,10,12-15"] [-e "O,H"] [-p PAD] [-v VAC] [-m]
```

| Option | Description |
|--------|-------------|
| `-i FILE` | Input `.xyz` file (required) |
| `-c "9 numbers"` | Lattice vectors `ax ay az bx by bz cx cy cz`. **Must be quoted.** |
| `-o FILE` | Output file (default `POSCAR`) |
| `-t TYPE` | `Cartesian` (default) or `Direct` |
| `-f RANGES` | 1-indexed atoms to freeze, e.g. `"1-5,10,12-15"` |
| `-e ELEMS` | Comma-separated elements to freeze, e.g. `"Au,O"` (union with `-f`) |
| `-v VAC` | Vacuum in angstrom on each side; box = extent + 2*VAC; molecule centred (implies `-m`) |
| `-m` | Centre the molecule in the cell (bounding-box midpoint placed at `(a+b+c)/2`) |
| `-p PAD` | Default box only (no `-v`): side = max coordinate + PAD (default `1.0`) |
| `-h` | Show help |

### Which cell is used?

| You pass | Cell |
|----------|------|
| `-c "..."` | your cell |
| nothing, extended-XYZ file with `Lattice="..."` | the cell from the file |
| `-v 10` | extent + 20 angstrom per axis, molecule centred |
| nothing else | max coordinate + `PAD` (default 1 angstrom) per axis |

If a cell is supplied or detected, `-v` is not used for the box size; the molecule is only centred.

## Example

Freeze the four Au atoms of a small CO/Au4 model and put it in a box with 10 angstrom of
vacuum on every side:

```bash
./xyz2poscar.sh -i examples/co_on_au4.xyz -v 10 -e "Au"
```

Messages printed to the terminal (stderr):

```text
Default box from vacuum: 22.8800 x 22.8800 x 23.1500 (extent + 2 x 10)
Molecule centred in the cell: all atoms translated by (11.440000, 11.440000, 10.000000)
Vacuum on each side of the molecule: x 10.0000, y 10.0000, z 10.0000
Wrote POSCAR (Cartesian coordinates, atoms grouped by element).
```

Resulting `POSCAR`:

```text
CO above a planar Au4 cluster (illustrative geometry)
1.0
       22.8800000000        0.0000000000        0.0000000000
        0.0000000000       22.8800000000        0.0000000000
        0.0000000000        0.0000000000       23.1500000000
Au C O 
4 1 1 
Selective dynamics
Cartesian
       10.0000000000       10.0000000000       10.0000000000  F F F  ! atom 1 (Au)
       12.8800000000       10.0000000000       10.0000000000  F F F  ! atom 2 (Au)
       10.0000000000       12.8800000000       10.0000000000  F F F  ! atom 3 (Au)
       12.8800000000       12.8800000000       10.0000000000  F F F  ! atom 4 (Au)
       11.4400000000       11.4400000000       12.0000000000  T T T  ! atom 5 (C)
       11.4400000000       11.4400000000       13.1500000000  T T T  ! atom 6 (O)
```

More patterns:

```bash
# your own cell, freeze atoms 1-4 and all oxygens, centre the molecule in that cell
./xyz2poscar.sh -i sample.xyz -c "17 0 0 0 15 0 0 0 35" -f "1-4" -e "O" -m

# extended-XYZ input: the cell is read from the Lattice="..." header
./xyz2poscar.sh -i extended.xyz -f "1-8"

# fractional coordinates, 12 angstrom vacuum
./xyz2poscar.sh -i sample.xyz -v 12 -t Direct -o POSCAR_direct
```

## Behaviour details

- **Shift to the positive octant** (only when `-m`/`-v` are not used): if any atom has a
  negative x, y or z, all atoms are translated by one common vector per axis so the minimum
  on that axis becomes 0. Axes without negative values are untouched. This is a rigid
  translation, so bond lengths and angles do not change.
- **Centring** (`-m`, `-v`): all atoms are translated so that the midpoint of the molecule's
  bounding box coincides with the cell centre. This replaces the positive-octant shift.
- **Default box** (no `-c`, no extended-XYZ lattice, no `-v`): orthorhombic, with each side
  equal to the largest coordinate along that axis (after the shift) plus `PAD`.
- **Atom order**: atoms are grouped by element in order of first appearance, keeping their
  relative order inside each group. Comments of the form `! atom N (El)` carry the original
  line number from the `.xyz` file. VASP ignores text after the flags.
- **Units**: coordinates and `-p`/`-v` are in the units of the input file (angstrom for VASP).
  The scaling factor is always `1.0`.

## Limitations

- Freezing is all-or-nothing (`F F F`); per-axis constraints require editing the `flag`
  assignment in the awk block.
- The automatic box is orthorhombic; use `-c` for other cell shapes.
- Centring uses the bounding-box midpoint, not the centre of mass, and acts in all three
  directions. For a periodic slab, supply your own cell with `-c`.
- Only the `Lattice="..."` entry of extended-XYZ headers is read; other properties are ignored.
- Atoms are not wrapped back into the cell, so check that a user-supplied cell is large enough.
- For MD runs choose the vacuum so the molecule stays far from its periodic images for the
  whole trajectory, not only at the starting geometry.

## Tests

```bash
bash tests/run_tests.sh
```

The tests cover freezing, positive-octant shift, default box, `-v`/`-m` centring, extended-XYZ
detection, Direct coordinates, geometry preservation and input validation. The same tests and
ShellCheck run on every push through GitHub Actions (`.github/workflows/ci.yml`).

## Repository layout

```text
xyz2poscar/
├── xyz2poscar.sh
├── README.md
├── LICENSE
├── .gitignore
├── examples/
│   ├── co_on_au4.xyz
│   ├── negative_coords.xyz
│   └── extended.xyz
├── tests/
│   └── run_tests.sh
└── .github/workflows/ci.yml
```

## License

MIT, see [LICENSE](LICENSE).

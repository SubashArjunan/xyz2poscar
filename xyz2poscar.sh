#!/usr/bin/env bash
#
# xyz2poscar.sh - Convert a .xyz file into a VASP POSCAR file with
#                 Selective Dynamics, allowing specific atoms to be frozen.
#
# Usage:
#   ./xyz2poscar.sh -i structure.xyz [-c "ax ay az bx by bz cx cy cz"] \
#                    [-o POSCAR] [-t Cartesian|Direct] \
#                    [-f "1-5,10,12-15"] [-e "O,H"] [-p PAD] [-v VAC] [-m]
#
# Options:
#   -i FILE     Input .xyz file (required)
#   -c "9 nums" Lattice vectors as 9 space-separated numbers:
#               ax ay az bx by bz cx cy cz
#               Optional. If omitted, the script tries to auto-detect a
#               cell from an extended-XYZ header; if none is found, it
#               builds a default box (see below).
#   -o FILE     Output file name (default: POSCAR)
#   -t TYPE     Coordinate type in output: Cartesian (default) or Direct
#   -f RANGES   1-indexed atom indices (order as in the .xyz file) to
#               freeze, e.g. "1-5,10,12-15". These atoms get F F F.
#   -e ELEMS    Comma-separated element symbols to freeze entirely,
#               e.g. "O,H". Combines with -f (union of both sets).
#   -v VAC      Vacuum (angstrom) on EACH side of the molecule. Builds an
#               orthorhombic box of (extent + 2*VAC) along x, y, z and
#               puts the molecule in the centre. Implies -m. Ignored for
#               the box size if -c or an extended-XYZ lattice is used
#               (the molecule is still centred in that cell).
#   -m          Centre the molecule in the cell (works with any cell:
#               -c, extended-XYZ lattice, or default box). The centre of
#               the molecule is the midpoint of its bounding box; it is
#               placed at the centre of the cell, (a+b+c)/2.
#   -p PAD      Used only for the default box when -v is NOT given:
#               box side = max coordinate (after shift) + PAD. Default 1.0
#   -h          Show this help.
#
# Automatic behaviour:
#   1) Without -m / -v: if any atom has a negative x, y or z coordinate,
#      ALL atoms are translated by one common vector per axis so that the
#      minimum coordinate on that axis becomes 0 (rigid translation,
#      geometry unchanged). Axes with no negative coordinates are left
#      untouched.
#   2) Cell selection, in order of priority:
#        a) -c given explicitly
#        b) extended-XYZ header on line 2, e.g.
#             Lattice="10.0 0.0 0.0 0.0 10.0 0.0 0.0 0.0 10.0" Properties=...
#        c) default orthorhombic box:
#             with -v : side = extent + 2*VAC
#             otherwise: side = max(x|y|z) + PAD  (max taken after shift)
#   3) With -m / -v the molecule is centred (rigid translation), which
#      replaces the shift in step 1.
#
# Notes:
#   * VASP requires atoms of the same species to be contiguous in POSCAR.
#     This script groups atoms by element (stable, first-seen-element
#     order), so the printed "atom N" comment refers to the ORIGINAL
#     line order in the .xyz file, which is what -f indices refer to.
#   * "Frozen" means fixed in all three directions (F F F). Edit the
#     flag assignment in the awk block below if you need per-axis control.
#   * Centring is done in all three directions. For a periodic slab you
#     usually want to supply your own cell with -c instead.
#   * Quote multi-word arguments: -c "17 0 0 0 15 0 0 0 35"

set -euo pipefail

OUTFILE="POSCAR"
COORD="Cartesian"
FREEZE_IDX=""
FREEZE_ELEM=""
XYZFILE=""
CELL=""
PAD="1.0"
VAC=""
CENTER=0

usage() {
    grep '^#' "$0" | sed -e '1d' -e 's/^# \{0,1\}//'
}

while getopts ":i:c:o:t:f:e:p:v:mh" opt; do
    case "$opt" in
        i) XYZFILE="$OPTARG" ;;
        c) CELL="$OPTARG" ;;
        o) OUTFILE="$OPTARG" ;;
        t) COORD="$OPTARG" ;;
        f) FREEZE_IDX="$OPTARG" ;;
        e) FREEZE_ELEM="$OPTARG" ;;
        p) PAD="$OPTARG" ;;
        v) VAC="$OPTARG"; CENTER=1 ;;
        m) CENTER=1 ;;
        h) usage; exit 0 ;;
        \?) echo "Unknown option: -$OPTARG" >&2; usage; exit 1 ;;
        :) echo "Option -$OPTARG requires an argument" >&2; exit 1 ;;
    esac
done

if [[ -z "$XYZFILE" ]]; then
    echo "Error: -i is required." >&2
    usage
    exit 1
fi

if [[ ! -f "$XYZFILE" ]]; then
    echo "Error: input file '$XYZFILE' not found." >&2
    exit 1
fi

# If no explicit cell was given, try to auto-detect an extended-XYZ
# Lattice="..." entry on the comment line (line 2 of the file).
if [[ -z "$CELL" ]]; then
    LATTICE_LINE=$(sed -n '2p' "$XYZFILE")
    if [[ "$LATTICE_LINE" =~ Lattice=\"([^\"]+)\" ]]; then
        CELL="${BASH_REMATCH[1]}"
        NCELL=$(wc -w <<< "$CELL")
        if [[ "$NCELL" -eq 9 ]]; then
            echo "Detected cell from extended-XYZ header: $CELL" >&2
        else
            echo "Warning: found a Lattice= tag but it doesn't have 9 numbers; ignoring it." >&2
            CELL=""
        fi
    fi
fi

if [[ -n "$CELL" ]]; then
    NCELL=$(wc -w <<< "$CELL")
    if [[ "$NCELL" -ne 9 ]]; then
        echo "Error: -c must contain exactly 9 numbers (three lattice vectors). Did you forget the quotes?" >&2
        exit 1
    fi
    if [[ -n "$VAC" ]]; then
        echo "Note: a cell was supplied/detected, so -v is not used for the box size; the molecule is only centred." >&2
    fi
fi

if [[ "$COORD" != "Cartesian" && "$COORD" != "Direct" ]]; then
    echo "Error: -t must be 'Cartesian' or 'Direct'." >&2
    exit 1
fi

# Expand a range spec like "1-5,10,12-15" into a flat, space-separated list
expand_ranges() {
    local spec="$1"
    local out=()
    local part start end i
    IFS=',' read -ra parts <<< "$spec"
    for part in "${parts[@]}"; do
        [[ -z "$part" ]] && continue
        if [[ "$part" == *-* ]]; then
            start="${part%-*}"
            end="${part#*-}"
            for ((i = start; i <= end; i++)); do out+=("$i"); done
        else
            out+=("$part")
        fi
    done
    echo "${out[*]}"
}

FREEZE_LIST=""
if [[ -n "$FREEZE_IDX" ]]; then
    FREEZE_LIST=$(expand_ranges "$FREEZE_IDX")
fi

awk -v cell="$CELL" -v coord="$COORD" -v pad="$PAD" -v vac="$VAC" -v center="$CENTER" \
    -v freeze_list="$FREEZE_LIST" -v freeze_elem="$FREEZE_ELEM" '
function build_inverse() {
    det = M[1]*(M[5]*M[9]-M[6]*M[8]) - M[2]*(M[4]*M[9]-M[6]*M[7]) + M[3]*(M[4]*M[8]-M[5]*M[7])
    if (coord == "Direct") {
        if (det == 0) {
            print "Error: lattice matrix is singular, cannot build Direct coordinates." > "/dev/stderr"
            exit 1
        }
        inv[1]=(M[5]*M[9]-M[6]*M[8])/det; inv[4]=(M[3]*M[8]-M[2]*M[9])/det; inv[7]=(M[2]*M[6]-M[3]*M[5])/det
        inv[2]=(M[6]*M[7]-M[4]*M[9])/det; inv[5]=(M[1]*M[9]-M[3]*M[7])/det; inv[8]=(M[3]*M[4]-M[1]*M[6])/det
        inv[3]=(M[4]*M[8]-M[5]*M[7])/det; inv[6]=(M[2]*M[7]-M[1]*M[8])/det; inv[9]=(M[1]*M[5]-M[2]*M[4])/det
    }
}
BEGIN {
    have_cell = 0
    if (cell != "") {
        split(cell, c, " ")
        for (i = 1; i <= 9; i++) M[i] = c[i] + 0
        have_cell = 1
    }
    n = split(freeze_list, farr, " ")
    for (i = 1; i <= n; i++) freeze_idx[farr[i]] = 1
    m = split(freeze_elem, earr, ",")
    for (i = 1; i <= m; i++) freeze_el[earr[i]] = 1
}
NR == 1 { natoms = $1; next }
NR == 2 { comment = $0; next }
NR > 2 && NF >= 4 {
    idx = NR - 2
    elem = $1
    order[idx]  = elem
    ax[idx] = $2 + 0; ay[idx] = $3 + 0; az[idx] = $4 + 0
    if (idx == 1) {
        minx = maxx = ax[idx]; miny = maxy = ay[idx]; minz = maxz = az[idx]
    } else {
        if (ax[idx] < minx) minx = ax[idx]; if (ax[idx] > maxx) maxx = ax[idx]
        if (ay[idx] < miny) miny = ay[idx]; if (ay[idx] > maxy) maxy = ay[idx]
        if (az[idx] < minz) minz = az[idx]; if (az[idx] > maxz) maxz = az[idx]
    }
    frozen[idx] = ((idx in freeze_idx) || (elem in freeze_el)) ? 1 : 0
    if (!(elem in seen)) { seen[elem] = 1; elemlist[++necat] = elem }
    count[elem]++
}
END {
    # extent of the molecule along each axis (translation invariant)
    ex = maxx - minx; ey = maxy - miny; ez = maxz - minz

    # shift that would bring the minimum coordinate to 0 on axes with negatives
    sx = (minx < 0) ? -minx : 0
    sy = (miny < 0) ? -miny : 0
    sz = (minz < 0) ? -minz : 0

    # --- choose the cell -------------------------------------------------
    if (!have_cell) {
        if (vac != "") {
            v = vac + 0
            M[1]=ex+2*v; M[2]=0;      M[3]=0
            M[4]=0;      M[5]=ey+2*v; M[6]=0
            M[7]=0;      M[8]=0;      M[9]=ez+2*v
            printf "Default box from vacuum: %.4f x %.4f x %.4f (extent + 2 x %s)\n", M[1], M[5], M[9], vac > "/dev/stderr"
        } else {
            p = pad + 0
            M[1]=maxx+sx+p; M[2]=0;         M[3]=0
            M[4]=0;         M[5]=maxy+sy+p; M[6]=0
            M[7]=0;         M[8]=0;         M[9]=maxz+sz+p
            printf "No cell given or detected: default orthorhombic box %.4f x %.4f x %.4f (max coordinate + %s)\n", M[1], M[5], M[9], pad > "/dev/stderr"
        }
    }
    build_inverse()

    # --- move the molecule -----------------------------------------------
    if (center + 0 == 1) {
        # put the bounding-box midpoint at the cell centre (a+b+c)/2
        cx0 = (M[1] + M[4] + M[7]) / 2
        cy0 = (M[2] + M[5] + M[8]) / 2
        cz0 = (M[3] + M[6] + M[9]) / 2
        tx = cx0 - (minx + maxx) / 2
        ty = cy0 - (miny + maxy) / 2
        tz = cz0 - (minz + maxz) / 2
        for (idx = 1; idx <= natoms; idx++) {
            ax[idx] += tx; ay[idx] += ty; az[idx] += tz
        }
        printf "Molecule centred in the cell: all atoms translated by (%.6f, %.6f, %.6f)\n", tx, ty, tz > "/dev/stderr"
        if (M[2]==0 && M[3]==0 && M[4]==0 && M[6]==0 && M[7]==0 && M[8]==0) {
            printf "Vacuum on each side of the molecule: x %.4f, y %.4f, z %.4f\n", (M[1]-ex)/2, (M[5]-ey)/2, (M[9]-ez)/2 > "/dev/stderr"
        }
    } else if (sx > 0 || sy > 0 || sz > 0) {
        for (idx = 1; idx <= natoms; idx++) {
            ax[idx] += sx; ay[idx] += sy; az[idx] += sz
        }
        printf "Negative coordinates found: all atoms shifted by (%.6f, %.6f, %.6f)\n", sx, sy, sz > "/dev/stderr"
    }

    # --- write POSCAR ----------------------------------------------------
    if (comment == "") comment = "Converted from XYZ by xyz2poscar.sh"
    print comment
    print "1.0"
    printf "%20.10f%20.10f%20.10f\n", M[1], M[2], M[3]
    printf "%20.10f%20.10f%20.10f\n", M[4], M[5], M[6]
    printf "%20.10f%20.10f%20.10f\n", M[7], M[8], M[9]

    line = ""
    for (i = 1; i <= necat; i++) line = line elemlist[i] " "
    print line
    line = ""
    for (i = 1; i <= necat; i++) line = line count[elemlist[i]] " "
    print line

    print "Selective dynamics"
    print coord

    for (i = 1; i <= necat; i++) {
        e = elemlist[i]
        for (idx = 1; idx <= natoms; idx++) {
            if (order[idx] != e) continue
            flag = frozen[idx] ? "F F F" : "T T T"
            if (coord == "Cartesian") {
                px = ax[idx]; py = ay[idx]; pz = az[idx]
            } else {
                px = ax[idx]*inv[1] + ay[idx]*inv[2] + az[idx]*inv[3]
                py = ax[idx]*inv[4] + ay[idx]*inv[5] + az[idx]*inv[6]
                pz = ax[idx]*inv[7] + ay[idx]*inv[8] + az[idx]*inv[9]
            }
            printf "%20.10f%20.10f%20.10f  %s  ! atom %d (%s)\n", px, py, pz, flag, idx, e
        }
    }
}
' "$XYZFILE" > "$OUTFILE"

echo "Wrote $OUTFILE ($COORD coordinates, atoms grouped by element)."

#!/usr/bin/env bash
# Minimal regression tests for xyz2poscar.sh
# Run from anywhere:  bash tests/run_tests.sh
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCRIPT="$ROOT/xyz2poscar.sh"
EX="$ROOT/examples"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

pass=0
fail=0

check() {
    local name="$1" expected="$2" actual="$3"
    if [[ "$expected" == "$actual" ]]; then
        pass=$((pass + 1))
        echo "ok   - $name"
    else
        fail=$((fail + 1))
        echo "FAIL - $name"
        echo "       expected: $expected"
        echo "       actual:   $actual"
    fi
}

# Diagonal of the lattice matrix of a POSCAR
box() {
    awk 'NR==3{a=$1} NR==4{b=$2} NR==5{c=$3} END{printf "%.4f %.4f %.4f", a, b, c}' "$1"
}

# Minimum x, y, z over all atoms of a POSCAR (atom lines start at line 10)
mins() {
    awk 'NR>=10{ if(NR==10||$1<x)x=$1; if(NR==10||$2<y)y=$2; if(NR==10||$3<z)z=$3 }
         END{printf "%.4f %.4f %.4f", x, y, z}' "$1"
}

NEG="$EX/negative_coords.xyz"

# 1. freezing by index and by element (union of both)
bash "$SCRIPT" -i "$NEG" -f 3 -e Au -o "$TMP/p1" 2>/dev/null
check "freeze by index + element" "3" "$(grep -c 'F F F' "$TMP/p1")"
check "remaining atoms are free"  "1" "$(grep -c 'T T T' "$TMP/p1")"

# 2. negative coordinates are shifted into the positive octant
bash "$SCRIPT" -i "$NEG" -o "$TMP/p2" 2>/dev/null
check "shift to positive octant" "0.0000 0.0000 0.0000" "$(mins "$TMP/p2")"

# 3. default box = max coordinate + 1
check "default box is max + 1" "5.0000 5.0000 6.5000" "$(box "$TMP/p2")"

# 4. -v : vacuum on each side, molecule centred
bash "$SCRIPT" -i "$NEG" -v 10 -o "$TMP/p4" 2>/dev/null
check "-v 10 box size"          "24.0000 24.0000 25.5000" "$(box "$TMP/p4")"
check "-v 10 vacuum below mol." "10.0000 10.0000 10.0000" "$(mins "$TMP/p4")"

# 5. explicit cell + centring
bash "$SCRIPT" -i "$NEG" -c "17 0 0 0 15 0 0 0 35" -m -o "$TMP/p5" 2>/dev/null
check "explicit cell is used" "17.0000 15.0000 35.0000" "$(box "$TMP/p5")"
check "centred in explicit cell" "6.5000 8.5000 17.2500" \
    "$(awk '/atom 1 /{printf "%.4f %.4f %.4f", $1, $2, $3}' "$TMP/p5")"

# 6. extended-XYZ lattice is auto-detected
bash "$SCRIPT" -i "$EX/extended.xyz" -o "$TMP/p6" 2>/dev/null
check "extended-XYZ lattice" "8.0000 8.0000 8.0000" "$(box "$TMP/p6")"

# 7. Direct coordinates fall inside [0,1)
bash "$SCRIPT" -i "$NEG" -v 10 -t Direct -o "$TMP/p7" 2>/dev/null
check "Direct coordinates in [0,1)" "0" \
    "$(awk 'NR>=10{for(i=1;i<=3;i++) if($i<0||$i>=1) n++} END{print n+0}' "$TMP/p7")"

# 8. rigid translation: distance between atoms 1 and 4 is preserved
d_in=$(awk 'NR==3{x1=$2;y1=$3;z1=$4} NR==6{x2=$2;y2=$3;z2=$4}
            END{printf "%.6f", sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)}' "$NEG")
d_out=$(awk '/atom 1 /{x1=$1;y1=$2;z1=$3} /atom 4 /{x2=$1;y2=$2;z2=$3}
             END{printf "%.6f", sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)}' "$TMP/p4")
check "geometry preserved (d[1-4])" "$d_in" "$d_out"

# 9. a malformed -c must be rejected
if bash "$SCRIPT" -i "$NEG" -c "1 2 3" -o "$TMP/p9" >/dev/null 2>&1; then
    check "malformed -c is rejected" "rejected" "accepted"
else
    check "malformed -c is rejected" "rejected" "rejected"
fi

echo
echo "$pass passed, $fail failed"
[[ "$fail" -eq 0 ]]

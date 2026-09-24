#!/bin/bash
# ============================================================================
# Which docked poses can actually be built in CHARMM?
#
#   ./check_poses.sh <system> <docked.pdbqt>
#
# AutoDock rotates the glycosidic torsions without sugar-specific torsional
# preferences, so some poses come back eclipsed: two ring hydrogens on adjacent
# sugars end up closer than 1.2 A once the hydrogens CHARMM needs are built.
# 02_build_glycan.py refuses those, and it is right to -- but the refusal is
# about ONE pose, not about the docking run, and the next pose down is usually
# fine and within a fraction of a kcal/mol.
#
# This tries every pose in the file and prints which are clean, so the choice is
# made on evidence instead of by assuming pose 1 is usable.
#
# It leaves build/<system>/ligand_af3.pdb set to the pose you pick, NOT to the
# last one tried -- rerun new_system.sh (or pose_from_docking.py) with the pose
# you want afterwards.
# ============================================================================
set -u
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"
SYS=${1:-}; POSEFILE=${2:-}
[ -n "${SYS}" ] && [ -n "${POSEFILE}" ] || { sed -n '2,20p' "$0"; exit 1; }
[ -f "systems/${SYS}.conf" ] || { echo "ERROR: no systems/${SYS}.conf"; exit 1; }
[ -f "${POSEFILE}" ] || { echo "ERROR: no such pose file: ${POSEFILE}"; exit 1; }
# shellcheck source=/dev/null
source "systems/${SYS}.conf"
LIG=${LIGAND:-X6}

NPOSE=$(grep -c '^MODEL' "${POSEFILE}"); [ "${NPOSE}" -eq 0 ] && NPOSE=1
SAVE=$(mktemp); cp "build/${SYS}/ligand_af3.pdb" "${SAVE}" 2>/dev/null || true

printf "%-6s %-12s %s\n" "pose" "affinity" "CHARMM build"
ok=""
for n in $(seq 1 "${NPOSE}"); do
  info=$(python3 pose_from_docking.py "${POSEFILE}" "build/${SYS}/ligand_af3.pdb" --pose "${n}" 2>&1) || {
    printf "%-6s %-12s %s\n" "${n}" "-" "pose extraction failed"; continue; }
  aff=$(echo "${info}" | sed -n 's/.*affinity \([-0-9.]*\).*/\1/p')
  out=$(LIGAND="${LIG}" python3 02_build_glycan.py "${SYS}" 2>&1)
  if echo "${out}" | grep -q "REFUSED"; then
    printf "%-6s %-12s %s\n" "${n}" "${aff:--}" \
      "REFUSED: $(echo "${out}" | sed -n 's/.*REFUSED: //p' | head -1)"
  else
    printf "%-6s %-12s %s\n" "${n}" "${aff:--}" \
      "ok   (closest H contact $(echo "${out}" | sed -n 's/.*closest non-bonded contact //p' | head -1))"
    [ -z "${ok}" ] && ok="${n}"
  fi
done

[ -s "${SAVE}" ] && cp "${SAVE}" "build/${SYS}/ligand_af3.pdb"; rm -f "${SAVE}"
echo ""
if [ -n "${ok}" ]; then
  echo "Best-affinity pose that builds: ${ok}"
  echo "Use it:  python3 pose_from_docking.py ${POSEFILE} build/${SYS}/ligand_af3.pdb --pose ${ok}"
  echo "         LIGAND=${LIG} python3 02_build_glycan.py ${SYS}"
else
  echo "NO pose builds. That is not a pose problem any more -- check the ligand"
  echo "topology in ligands/${LIG}_charmm/ matches this molecule."
fi

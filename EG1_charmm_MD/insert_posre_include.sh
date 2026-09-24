#!/bin/bash
# Inserts a `#ifdef <DEFINE>\n#include "<itp>"\n#endif` block into a GROMACS
# .top file, right after the [ atoms ] section of the protein moleculetype --
# found automatically as the moleculetype with the most atoms (water/ion
# moleculetypes have only 1-4 atoms, the protein has hundreds+), rather than
# assuming it's literally named "Protein" (acpype's naming may differ from
# pdb2gmx's).
#
# Usage: ./insert_posre_include.sh system.top POSRES posre.itp
#        ./insert_posre_include.sh system.top POSRES_BB posre_backbone.itp
#
# Edits the .top in place (after writing a .bak) and prints the inserted
# region back so you can visually confirm placement.
set -e

TOP=$1
DEFINE=$2
ITP=$3

if [ -z "${TOP}" ] || [ -z "${DEFINE}" ] || [ -z "${ITP}" ]; then
  echo "Usage: $0 <system.top> <DEFINE_NAME> <itp_file>"
  exit 1
fi
if [ ! -f "${TOP}" ]; then
  echo "ERROR: ${TOP} not found"
  exit 1
fi

total_lines=$(wc -l < "${TOP}")

# Line numbers of every "[ atoms ]" header, and every section header of any
# kind. Using a plain while-read loop instead of `mapfile` (bash 4+ only) --
# this way it also works on macOS's stock bash 3.2, not just Linux/HPC bash.
ATOM_LINES=()
while IFS= read -r ln; do ATOM_LINES+=("${ln}"); done \
  < <(grep -n "^\[[ \t]*atoms[ \t]*\]" "${TOP}" | cut -d: -f1)
ALL_SECTIONS=()
while IFS= read -r ln; do ALL_SECTIONS+=("${ln}"); done \
  < <(grep -n "^\[" "${TOP}" | cut -d: -f1)

if [ "${#ATOM_LINES[@]}" -eq 0 ]; then
  echo "ERROR: no [ atoms ] sections found in ${TOP}"
  exit 1
fi

best_header=-1
best_count=-1
best_end=-1

for atoms_header in "${ATOM_LINES[@]}"; do
  # atoms_end = next section header after this one, or one past EOF if none
  atoms_end=$((total_lines + 1))
  for s in "${ALL_SECTIONS[@]}"; do
    if [ "${s}" -gt "${atoms_header}" ]; then
      atoms_end=${s}
      break
    fi
  done
  # count non-blank, non-comment lines strictly between the header and atoms_end
  count=$(sed -n "$((atoms_header + 1)),$((atoms_end - 1))p" "${TOP}" | grep -v '^[[:space:]]*;' | grep -c '[^[:space:]]' || true)
  if [ "${count}" -gt "${best_count}" ]; then
    best_count=${count}
    best_header=${atoms_header}
    best_end=${atoms_end}
  fi
done

echo "[ atoms ] sections found (line, atom count):"
for atoms_header in "${ATOM_LINES[@]}"; do
  atoms_end=$((total_lines + 1))
  for s in "${ALL_SECTIONS[@]}"; do
    if [ "${s}" -gt "${atoms_header}" ]; then
      atoms_end=${s}
      break
    fi
  done
  count=$(sed -n "$((atoms_header + 1)),$((atoms_end - 1))p" "${TOP}" | grep -v '^[[:space:]]*;' | grep -c '[^[:space:]]' || true)
  echo "  line ${atoms_header}: ${count} atoms"
done
echo "Inserting after line $((best_end - 1)) (${best_count} atoms -- the largest, assumed to be the protein)."
echo "If that's wrong (e.g. multiple protein copies similar in size to something else), stop and check ${TOP} manually."

insert_at=$((best_end - 1))
cp "${TOP}" "${TOP}.bak"
{
  sed -n "1,${insert_at}p" "${TOP}"
  echo "#ifdef ${DEFINE}"
  echo "#include \"${ITP}\""
  echo "#endif"
  sed -n "$((insert_at + 1)),\$p" "${TOP}"
} > "${TOP}.new"
mv "${TOP}.new" "${TOP}"

echo ""
echo "Wrote ${TOP} (backup saved as ${TOP}.bak). Inserted region:"
sed -n "$((insert_at - 1)),$((insert_at + 4))p" "${TOP}"

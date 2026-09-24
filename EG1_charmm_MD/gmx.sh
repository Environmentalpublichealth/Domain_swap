# ============================================================================
# Resolve a working `gmx` and expose it as ${GMX}. Sourced, not executed:
#
#     source "${SCRIPT_DIR}/gmx.sh"
#     ${GMX} grompp -f ... ; ${GMX} mdrun ...
#
# On this cluster GROMACS lives in an Apptainer container and `gmx` is NOT on
# the container's PATH -- the arch-specific full path has to be given, or you
# get the misleading
#     FATAL: "gmx": executable file not found in $PATH
# (lab wiki: Domain_swap/wiki/GROMACS).
#
# --nv passes the host's NVIDIA driver into the container. Without it the
# container starts fine and then reports no GPU, which looks like a cluster
# problem but is a missing flag.
#
# If GMX is already set in the environment it is left alone, so a different
# cluster (e.g. a plain `module load gromacs`) needs no edit here:
#     GMX="gmx_mpi" ./prep_system.sh ...
# ============================================================================
if [ -z "${GMX:-}" ]; then
  if [ -n "${CONTAINER:-}" ] && [ -e "${CONTAINER}" ]; then
    command -v apptainer >/dev/null 2>&1 && APPT=apptainer || APPT=singularity
    command -v "${APPT}" >/dev/null 2>&1 || {
      echo "FATAL: neither apptainer nor singularity is on PATH."; exit 1; }
    # Bind the project directory so the container can see the working files.
    # Apptainer mounts $HOME by default but NOT /engrfs.
    # Same form as the lab's working GPU job: --nv for the driver, an explicit
    # host:container bind of the project tree (Apptainer mounts $HOME but not
    # /engrfs, so without this the container cannot see the run directory).
    BIND_ROOT="${BIND_ROOT:-$(cd "${PIPELINE_ROOT:-$PWD}" && pwd)}"
    GMX="${APPT} exec --nv --bind ${BIND_ROOT}:${BIND_ROOT} ${CONTAINER} ${GMX_IN_CONTAINER}"
  elif command -v gmx >/dev/null 2>&1; then
    GMX="gmx"                       # thread-MPI build on PATH
  elif command -v gmx_mpi >/dev/null 2>&1; then
    GMX="gmx_mpi"
  else
    echo "FATAL: no GROMACS found."
    echo "  container looked for : ${CONTAINER:-<unset>}"
    echo "  and no gmx/gmx_mpi on PATH."
    echo "Set CONTAINER in config.sh, or GMX=... in your environment."
    exit 1
  fi
fi

# Prove it runs before a job spends an hour finding out it does not. The version
# banner also says whether this build has GPU support at all.
gmx_check () {
  local v
  v=$(${GMX} --version 2>&1) || { echo "FATAL: '${GMX}' does not run:"; echo "${v}" | tail -5; return 1; }
  echo "gmx        : ${GMX}"
  echo "version    : $(echo "${v}" | grep -m1 -i 'GROMACS version' | sed 's/^ *//')"
  echo "GPU support: $(echo "${v}" | grep -m1 -i 'GPU support' | sed 's/.*: *//')"
  return 0
}

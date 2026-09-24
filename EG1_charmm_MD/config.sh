# ============================================================================
# EG1 domain-swap complexes -- CHARMM36m MD, one place for every setting.
#
# Sourced by every script here. Edit THIS file, never a copy inside a script:
# two files that disagree about the protocol is how two systems quietly end up
# incomparable.
#
# PART 1 is the cluster (WashU engineering RCC). PART 2 is the simulation
# protocol, copied unchanged from XYN_charmm_MD so the EG1 runs can be compared
# with the XYN runs directly.
# ============================================================================

# ============================================================ PART 1: CLUSTER
# Taken from a working job on this cluster (Docking/run_GNINA_4HK8_validation.slurm).
# If you submit under a different PI, change it; check yours with:
#   sacctmgr -n show assoc user=$USER format=account%30
ACCOUNT="${ACCOUNT:-engr-lab-joshua.yuan}"

# --- GROMACS comes from a container on this cluster, not a module.
# The container holds several builds, one per CPU instruction set, and `gmx`
# alone is NOT on the container's PATH -- the full path must be given.
# See the lab wiki: Domain_swap/wiki/GROMACS
CONTAINER="${CONTAINER:-/engrfs/project/joshua.yuan/Haina/software/containers/gromacs.sif}"
GMX_ARCH="${GMX_ARCH:-avx_512}"       # avx_512 | avx2_256 | avx_256 | sse4.1
                                      # check yours: lscpu | grep -o 'avx512[a-z]*' | head -1
GMX_IN_CONTAINER="${GMX_IN_CONTAINER:-/usr/local/gromacs/${GMX_ARCH}/bin/gmx}"

# --- Partitions, from the lab's own working GPU job.
PARTITION_CPU="${PARTITION_CPU:-general-cpu}"
PARTITION_GPU="${PARTITION_GPU:-general-gpu}"

# --- GPU request. NOT the usual "gpu:1": this cluster hands out SHARDS, i.e.
# fractions of a physical GPU, and the working example asks for
#     --gres="shard:a6000:1"        (the lab's GNINA example)
# We ask for an A100 SXM4 instead:
#     --gres="shard:a100-sxm4:1"
# A plain --gres=gpu:1 may be rejected or queue forever if the partition only
# advertises shards. See what is actually offered with:
#     sinfo -p general-gpu -o "%P %G %l"
# A shard is enough for one GROMACS replicate; ask for more only if a run is
# demonstrably GPU-bound.
GRES="${GRES:-shard:a100-sxm4:1}"
MEM="${MEM:-16G}"                     # the GNINA example used 8G; MD writes more

# --- Job shapes. Production runs on the GPU; prep is CPU-only and short.
GPU_CPUS="${GPU_CPUS:-8}"             # OpenMP threads beside the GPU
PREP_CPUS="${PREP_CPUS:-4}"
MAXCHAIN="${MAXCHAIN:-12}"            # self-resubmits before giving up

# --- Walltime: NOT REQUESTED. No --time is passed, so each job gets the
# partition's default limit.
#
# run_md.batch does not assume that limit. At startup it asks the scheduler what
# this job was actually granted (scontrol show job) and sets mdrun -maxh half an
# hour below it, so mdrun stops cleanly and checkpoints while the script is still
# alive to submit the next link. That matters: if mdrun is still running when
# SLURM SIGKILLs the job, the script dies too, the chain stops silently, and you
# find a replicate that just stopped with no error anywhere.
#
# MAXH is only the fallback for when the scheduler will not say (an UNLIMITED or
# unreadable limit). Being wrong here is cheap: -cpt 15 means at most 15 minutes
# of lost work, and the next link resumes from the checkpoint.
MAXH="${MAXH:-23.5}"
WALL_H=""                             # empty on purpose -- see above

# --- mdrun GPU offload.
# -nb/-pme/-bonded gpu is the safe, always-supported set. `-update gpu` is a
# further speedup but is refused in some builds/ensembles, so it is OFF by
# default; set GPU_UPDATE=1 to try it and read what mdrun says at startup.
GPU_FLAGS="${GPU_FLAGS:--nb gpu -pme gpu -bonded gpu}"
GPU_UPDATE="${GPU_UPDATE:-0}"
[ "${GPU_UPDATE}" = "1" ] && GPU_FLAGS="${GPU_FLAGS} -update gpu"

# =========================================================== PART 2: PROTOCOL
# Identical to XYN_charmm_MD/config.sh. Do not drift from it without deciding
# that the two projects are no longer being compared.
WATER_MODEL=tip3p          # CHARMM-modified TIP3P in this port
SALT_CATION=SOD            # CHARMM ion names. NOT Amber's NA/CL.
SALT_ANION=CLA
BOX_TYPE=dodecahedron
BOX_PAD=1.2                # nm, solute-to-box-edge
TEMP_K=300
SALT_MM=150                # bulk NaCl on top of neutralisation

NVT_PS=500                 # heavy atoms restrained (+ ligand, in holo)
NPT_NS=2                   # backbone restrained; ligand free
PROD_NS=200                # unrestrained production
FRAME_PS=10                # -> 20001 frames per replicate

N_REPS=2                   # independent replicates, reported individually

# --- Systems.
# Leave this EMPTY and every script treats "all systems" as "every
# systems/*.conf", so ./new_system.sh is the only place a system is declared and
# there is no second list to forget to update.
# Set it to a space-separated subset only to deliberately narrow what "all" means.
SYSTEMS=""

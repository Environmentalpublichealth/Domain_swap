# AI-HDX step 1: HHblits profiles for the 1XYN ↔ 1KS5 swap mutants

This folder makes the HHblits profiles (`.hhm` files) that AI-HDX needs as input. The protocol follows the AI-HDX documentation [`Documentations/MSA_embedding.md`](https://github.com/Environmentalpublichealth/AI-HDX/blob/main/Documentations/MSA_embedding.md): one protein per FASTA file, `hhblits` with default settings against UniRef30 2020_06, and each result saved as `hhm_data/<name>.hhm`. Everything runs as SLURM jobs on the WashU cluster.

**Sequences (12):**
- the 10 region-swap mutants: M02, M03, M05, M08 and M13, each in both directions:
  - `1KS5_Mxx` = a 1XYN segment put into 1KS5
  - `1XYN_Mxx` = a 1KS5 segment put into 1XYN
- the 2 wild types, `1XYN_WT` and `1KS5_WT`, as references.

---

## What's in the folder

| Path | What it is |
|---|---|
| `input/1KS5_1XYN_region_swap_10_sequences.fasta` | The 10 swap sequences (source) |
| `input/WT_reference.fasta` | Wild-type 1XYN (178 aa) and 1KS5 (223 aa), taken from the PDB crystal structures |
| `fasta/<name>.fasta` | **Ready-made inputs**, one sequence per file (made by `make_inputs.py`) |
| `inputs.tsv` | List of the sequences with their lengths. Row *i* = array task *i* |
| `config.sh` | **Settings to check**: SLURM account and partition, conda environment name, database location |
| `make_inputs.py` | Rebuilds `fasta/` and `inputs.tsv` from `input/` |
| `submit.sh` | Submits the HHblits job array |
| `slurm/hhblits_array.sbatch` | One array task = one sequence |
| `slurm/extract_db.sbatch` | One-time extraction of the database archive, **only if needed** (step 2) |
| `slurm/common.sh` | Shared helpers (conda activation, database check) |
| `check_outputs.sh` | Checks that every `.hhm` exists and has the right length |

Created when you run:
- `hhm_data/`: **the results** (`<name>.hhm`)
- `hhr/`: HHblits search reports (`<name>.hhr`); useful for troubleshooting, not needed by AI-HDX
- `logs/`: SLURM logs

---

## Step 0: Copy the folder to the cluster

Copy the whole `AI-HDX_hhblits/` folder into your own space on the cluster, for example your home or scratch folder. All commands below are run from inside that folder:

```bash
cd /path/to/AI-HDX_hhblits
```

## Step 1: Check the conda environment

```bash
conda env list                 # "hhblits" should be listed
conda activate hhblits
hhblits -h | head -3           # should print the HHblits version
conda deactivate
```

If the environment is missing, create it as in the AI-HDX documentation:
```bash
conda create -n hhblits
conda activate hhblits
conda install -c conda-forge -c bioconda hhsuite
```
If you use a different environment name, change `CONDA_ENV` in `config.sh`.

The job scripts activate the environment themselves. You don't need to activate it before submitting.

## Step 2: Make sure the database is extracted

The database is shared at:

```
/engrfs/project/joshua.yuan/Haina/software/database/UniRef30_2020_06_hhsuite.tar.gz
```

That is a **compressed archive**, and HHblits can't read it directly. It needs the extracted files: `UniRef30_2020_06_a3m.ffdata`, `_hhm.ffdata` and `_cs219.ffdata`, each with a matching `.ffindex` file.

**First, check whether someone has already extracted it:**

```bash
find /engrfs/project/joshua.yuan/Haina/software/database -maxdepth 3 -name '*_a3m.ffdata'
```

- **A file is listed:** nothing to do. With the default `config.sh`, the jobs find it automatically.
- **Nothing is listed:** the database has to be extracted once. **Don't extract it into Haina's folder without asking her.** Extract into your own space instead:
  1. In `config.sh`, set `DB_DIR` to a folder you own with plenty of free space. The extracted database is several times the size of the archive, and the script checks for at least 4× the archive size before starting. For example:
     ```bash
     DB_DIR=/path/to/your/space/databases/UniRef30_2020_06
     ```
  2. Extract, which can take several hours. The job runs on `general-cpu` under the lab account:
     ```bash
     mkdir -p logs
     sbatch slurm/extract_db.sbatch
     ```
  3. When it finishes, `logs/uniref30_extract_<jobid>.out` should list the `.ffdata` files.

Note: the AI-HDX documentation writes `tar -xvfz`. That flag order is wrong, because `f` must come right before the file name. The script uses `tar -xzf`.

## Step 3: Inputs (already prepared)

`fasta/` and `inputs.tsv` are already in the folder, so you can go straight to step 4. To add or change sequences:
1. Put FASTA files in `input/`. The name is the header text before the first `|` or space, for example `>1KS5_M02 | ...` becomes `1KS5_M02`. Names must be unique.
2. Rebuild:
   ```bash
   python3 make_inputs.py
   ```
It checks for empty sequences, non-standard letters and duplicate names, then rewrites `fasta/` and `inputs.tsv`.

## Step 4: Submit

```bash
bash submit.sh
```

This checks that the database is usable, then submits one array task per sequence (12 tasks). It uses the partition and account set in `config.sh`, which are the lab's usual settings:

| Setting in `config.sh` | Default | Notes |
|---|---|---|
| `PARTITION_CPU` | `general-cpu` | HHblits uses CPUs only, so it runs on the CPU partition |
| `ACCOUNT` | `engr-lab-joshua.yuan` | The lab's SLURM account. Check yours with `sacctmgr -n show assoc user=$USER format=account%30` |

Each task uses 1 node, 8 CPUs and 32 GB of memory. No time limit is requested, so the partition's default applies, as in the lab's other jobs on this cluster.

- **To change resources for one submission:** add options to `submit.sh`; they override the defaults. For example, `bash submit.sh --mem=64G`.
- **To change them permanently:** edit the `#SBATCH` lines at the top of `slurm/hhblits_array.sbatch`, or `config.sh` for the partition and account.

**Settings follow the AI-HDX documentation.** The command is:

```
hhblits -i fasta/<name>.fasta -ohhm hhm_data/<name>.hhm -o hhr/<name>.hhr -d <database> -cpu 8
```

This uses HHblits' default settings (2 search iterations, default E-value and filters), as AI-HDX does. `-cpu` only changes speed. **Don't add or change other options**: AI-HDX was trained on profiles made with these defaults.

## Step 5: Monitor and check

```bash
squeue -u $USER                    # running/pending jobs
tail logs/hhblits_<jobid>_<task>.out
bash check_outputs.sh              # when all tasks are finished
```

`check_outputs.sh` prints one line per sequence:

```
name         length   LENG sequences_passed_filter                  NEFF  status
1KS5_M02        223    223 1111 of 2026                              6.2  OK
```

- **LENG** must equal the sequence length. Each task already checks this before keeping the `.hhm`.
- **NEFF** is the effective number of distinct sequences in the alignment, a measure of alignment depth. Very low values (close to 1) mean few homologs were found, and AI-HDX predictions for that protein will be less reliable. Write down any sequence with a much lower NEFF than its wild type.

Any `MISSING` sequence: read `logs/hhblits_<jobid>_<task>.out`, fix the problem, and run `bash submit.sh` again. Finished sequences are skipped automatically.

## Step 6: Hand back the results

The results are the 12 files in `hhm_data/`. Copy that folder back to the shared drive, for example:

```bash
rsync -av <user>@<cluster>:/path/to/AI-HDX_hhblits/hhm_data .
```

These `.hhm` files are the MSA input for the AI-HDX embedding and prediction steps.

---

## Troubleshooting

| Message in the log | Fix |
|---|---|
| `FATAL: conda not on PATH` | Log in again, or run `conda init bash` once and start a new session |
| `FATAL: cannot activate conda env 'hhblits'` | Check the name with `conda env list`; update `CONDA_ENV` in `config.sh` |
| `FATAL: hhblits not in .../bin` | Install hhsuite into the environment (step 1) |
| `FATAL: no extracted database ... under DB_DIR=` | Step 2: point `DB_DIR` at the extracted copy, or extract it |
| `FATAL: ..._cs219.ffdata/.ffindex missing -- incomplete extraction?` | The extraction stopped part-way. Delete the partial files and run `extract_db.sbatch` again |
| `FATAL: <name>.hhm LENG=..., expected ...` | HHblits failed part-way; read the rest of the log. The partial result is kept as `hhm_data/<name>.hhm.tmp` |
| Job killed for memory (`OUT_OF_MEMORY`) | Resubmit with more memory: `bash submit.sh --mem=64G` |
| Job killed for time (`TIMEOUT`) | Resubmit with a longer limit: `bash submit.sh --time=24:00:00` |
| `Invalid account or account/partition combination` | Check your account with `sacctmgr -n show assoc user=$USER format=account%30` and set `ACCOUNT` in `config.sh` |

## Notes on the sequences

- Swap sequences are full length: each swap replaces a segment with one of the same length from the other protein. 1KS5 mutants are 223 aa and 1XYN mutants are 178 aa.
- The header ranges for **M03** in `input/1KS5_1XYN_region_swap_10_sequences.fasta` are wrong. They say `1KS5:44-51 ↔ 1XYN:23-30`; the sequences actually swap **1KS5:42-52 ↔ 1XYN:22-32**. The sequences are correct, and `make_inputs.py` writes only the name to each FASTA header, so the wrong ranges never reach the `.hhm` files.

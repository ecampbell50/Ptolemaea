# Ptolemaea Containerisation — Working Log

Goal: containerise the Ptolemaea pipeline with Apptainer so it runs reproducibly on
the QUB HPC (Kelvin), replacing Prokka with pyrodigal, and ultimately exposing a
single `ptolemaea <genomes_dir>` entry point.

## Key decisions
- **Apptainer**, not Docker (no root/daemon on HPC). Available on Kelvin via
  `module load apps/apptainer/1.3.4`.
- **One `.sif` per tool** (BioContainers / `apptainer pull`), NOT one conda env per
  tool in a single image — chosen because conda solves have been painful and Kelvin
  has no full `fakeroot` (so `apt-get` during a build is unreliable). Pulling
  pre-built images sidesteps both problems.
- Tools: pyrodigal (replaces Prokka), padloc 2.0.0, defense-finder, blast 2.15.0,
  hmmer 3.4. Final analysis (pandas) — TBD whether its own image or host python.
- Build/work in **scratch**, not home (home near-full). Cache redirected via
  `APPTAINER_CACHEDIR` / `APPTAINER_TMPDIR` in `~/.bashrc`.
- Keep this recipe/log in git so a scratch purge can't destroy progress — `.sif`s are
  rebuildable from the pull commands.
- End-user vision: Tier 1 = `ptolemaea genomes/` wrapper over the 5 `.sif`s (thesis
  target). Tier 2 = single self-contained image (future / publication).

## Environment facts (Kelvin)
- Apptainer 1.3.4 (module `apps/apptainer/1.3.4`; also 1.1.2 available)
- `fakeroot command not found` → builds use root-mapped namespace; `apt-get` installs
  in a build may fail. Pulling pre-built images is unaffected.
- `APPTAINER_CACHEDIR` and `APPTAINER_TMPDIR` set in `~/.bashrc` (dirs created).

## Plan (tick off as we go)
- [x] 0. Scratch working dir set up; cache env vars confirmed active
- [x] 1. Pull all 5 tool images as `.sif` (hmmer ✅ done as a test)
- [x] 2. Smoke-test each `.sif` (pyrodigal/blast/hmmer/defense-finder OK; padloc DB-write
       error is expected, handled in Step 4)
- [x] 3. Validate pyrodigal → PADLOC: PROVEN — padloc parsed pyrodigal faa/gff with NO shim,
       called real systems (e.g. PDC-S13). System set CONFIRMED identical to Prokka-based run.
- [x] 4. Databases done. PADLOC DB v2.0.0 in padloc_data. DefenseFinder runs with
       thesis-exact models (defense-finder-models 2.0.2 + CasFinder 3.1.0) copied from
       ~/.macsyfinder/models into scratch df_models. Output columns = standard 2.0.1 format.

## DESIGN REQUIREMENT: central, easy-to-change version pinning
All tool versions must be set in ONE obvious place (e.g. a `versions.sh`/config block or
vars at the top of the wrapper), so a user can bump a single value to upgrade. Current
pinned set (the paper/thesis versions — keep as default):
  PYRODIGAL=3.7.1  PADLOC=2.0.0 (DB v2.0.0)  BLAST=2.16.0  HMMER=3.4
  DEFENSEFINDER=2.0.1  +  DF_MODELS=2.0.2 / CasFinder 3.1.0
CRITICAL COUPLING: DefenseFinder image version and models version MUST move together
(DF 2.0.x <-> models 2.x; DF 3.0.0 <-> models 3.1.0). So "switch DF 2.0.1 -> 3.0.0" means
ALSO swapping the .sif tag AND re-running `defense-finder update` (latest models) — the
config/docs must make this pairing explicit so a user doesn't bump one and not the other.
- [x] 5. Rewrote `pipeline_functions.sh` + new `ptolemaea.config` (central version/path
       config). All tools via `apptainer exec`; run_prokka->run_pyrodigal (contig-rename
       dropped); padloc data overlay bind; DF --models-dir; blast via image; pandas steps
       via `ptol_python` (pandas.sif); add_genome_id_to_faa rewritten in awk (host).
       TODO before testing: (a) pull pandas.sif; (b) fix callers Ptolemaea.sh +
       defence_pipeline_array.sh (run_prokka->run_pyrodigal, 01_prokka->01_pyrodigal,
       drop module loads). PERF: every exec unpacks SIF (no squashfuse) -> consider
       enabling squashfuse or node-local APPTAINER_TMPDIR for the 700-genome SLURM run.
       NOTE: user re-wrote the file as `pipeline_functions_v2.sh` with their own teaching
       comments (CPUS default 8 to match sbatch); finished by assistant. When reviewed,
       `mv pipeline_functions_v2.sh pipeline_functions.sh` to make it canonical.
- LAYOUT: containerised work kept SEPARATE from the live pipeline. Original conda-based
  `scripts/pipeline_functions.sh` + `scripts/Ptolemaea.sh` are UNTOUCHED (still used for
  current runs). Containerised versions live in `singularity_scripts/`
  (pipeline_functions_v2.sh + ptolemaea.config). pandas.sif pulled (pandas 2.2.1).
- [x] 6. FULL END-TO-END SUCCESS on test genome via Ptolemaea_singularity.sh. All 5 stages
       ran in containers -> output/05_consensus/GCF_046031115.2_defenceprofile.csv with 65
       defence proteins (AGREE 12, RESOLVED 2, SINGLE 19, BLAST 17, MAPPING 15). pandas
       consensus step ran through pandas.sif via ptol_python. Harmless warnings: DF "not
       latest model" (we pinned 2.0.2 on purpose); BLAST "examine 5+ matches" (max_target_seqs 1).
       ~13 min wall for one genome (incl. SIF sandbox-conversion overhead each exec).
- 2026-06-11: VALIDATION CLOSED. diff of containerised vs original-thesis consensus for
  GCF_046031115.2: biologically IDENTICAL (same 65 proteins, same systems/statuses/outcomes).
  Only difference = protein ID scheme (pyrodigal real contig names NZ_CM170101.1_60 vs
  Prokka locus tags BFNECFPA_00061) + a cosmetic trailing space in the OLD prokka column.
  Pyrodigal IDs are arguably better (traceable to assembly). Containerised pipeline
  reproduces the paper.
- [ ] 7. Wire into SLURM array script (`apptainer exec` per task) + test a small batch
- [ ] 8. Write the `ptolemaea <genomes_dir>` wrapper front door
- [ ] 9. Update README (prereqs → Apptainer + images; Prokka → pyrodigal)

## IMPORTANT: bind-mount requirement
Kelvin does NOT auto-mount `/mnt/scratch2` into containers (only $HOME, /tmp). EVERY
`apptainer exec` that touches files on scratch MUST include `--bind /mnt/scratch2`,
else the tool reports "No such file or directory" for files that clearly exist.
This goes into every pipeline call and the final wrapper script.

PORTABILITY (Tier 2 / other systems): `/mnt/scratch2` is Kelvin-specific. Do NOT hardcode
it in anything meant to be portable. The wrapper script should bind the *working directory
the user passes in* (e.g. `--bind "$(realpath "$WORKDIR")"`) and/or honour the
`APPTAINER_BIND`/`APPTAINER_BINDPATH` env var, so each site can point it at their own
storage. Hardcoded scratch path is fine for our own runs only.

## Running log (problems & fixes)
- 2026-06-11: Confirmed `apptainer build --fakeroot test.sif docker://hello-world`
  works on Kelvin despite no fakeroot binary (root-mapped namespace fallback).
- 2026-06-11: `apptainer pull hmmer.sif docker://quay.io/biocontainers/hmmer:3.4--hdbdd923_1`
  succeeded; `hmmsearch -h` runs HMMER 3.4 from the .sif. `squashfuse not found` +
  "Converting SIF to temporary sandbox" warnings are harmless.
- 2026-06-11: Step 0 done. Cache vars confirmed -> /mnt/scratch2/users/40204129/apptainer/{cache,tmp}.
  Images dir: /mnt/scratch2/users/40204129/Ptolemaea_Singularity/images (hmmer.sif moved here).
- BioContainer tags are sourced live from the Quay API (build-hash suffixes change), e.g.:
  `curl -s "https://quay.io/api/v1/repository/biocontainers/<tool>/tag/?limit=15&onlyActiveTags=true" | python -c "import sys,json;[print(t['name']) for t in json.load(sys.stdin)['tags']]"`
- 2026-06-11: PINNED IMAGE TAGS (Step 1):
  - hmmer:         3.4--hdbdd923_1            (already pulled)
  - padloc:        2.0.0--hdfd78af_1
  - pyrodigal:     3.7.1--py312h247cb63_1     (py flavour irrelevant, it's a CLI)
  - defense-finder:2.0.1--pyhdfd78af_0        (NOT 3.0.0 newest — 3.x is a major bump that
                   may change genes.tsv columns/models and break create_defence_profile_direct.py;
                   2.0.1 matches the validated pipeline. Upgrade = deliberate revalidation later.)
  - blast:         2.16.0--h66d330f_5         (2.15.0 unavailable on registry; minor bump,
                   output stable for our fixed blastp/makeblastdb params)
- 2026-06-11: Step 2 smoke tests. Confirmed: pyrodigal v3.7.1, blast 2.16.0+, HMMER 3.4.
  defense-finder runs (version cmd is `defense-finder version`, not `--version`).
  padloc `--version` FAILS with: `mkdir: can't create directory '/usr/local/bin/../data':
  Read-only file system` — padloc defaults to storing its DB INSIDE the read-only image.
  Expected; fix in Step 4 = point padloc's data dir at a writable scratch bind-mount
  (PADLOC_DB / --data). Not a broken image.
- 2026-06-11: Step 3 — first pyrodigal run FAILED with "No such file or directory" for a
  file that existed. Cause: scratch not auto-mounted (see bind-mount note above). Fixed by
  adding `--bind /mnt/scratch2` to the exec call.
- 2026-06-11: pyrodigal output is PADLOC-compatible by inspection: faa first-token IDs
  (e.g. `NZ_CM170101.1_1`) MATCH gff `ID=` attributes. Uses real contig names (better than
  Prokka's renamed contig00001). 4850 proteins on the E. coli test genome. Residual risk:
  PADLOC must ignore the ` # ... ` metadata after the first token in faa headers — strip
  with sed only if it complains.
- 2026-06-11: PADLOC read-only data-dir error SOLVED. padloc tries to own `/usr/local/data`
  (= `/usr/local/bin/../data`) inside the read-only image and errexits on mkdir at startup
  (before parsing args, so even --help failed). FIX: overlay a writable scratch dir onto
  that path:
    mkdir -p .../padloc_data
    apptainer exec --bind /mnt/scratch2 \
      --bind .../padloc_data:/usr/local/data padloc.sif padloc ...
  padloc --help now runs. DB downloads into padloc_data on scratch (persists, survives
  outside the .sif). Example inputs live at /usr/local/test inside the image.
- 2026-06-11: padloc --db-update SUCCESS -> database **v2.0.0** downloaded/compiled into
  scratch padloc_data (matches the DB version of the original validated pipeline). BusyBox
  `rm`/`mkdir` warnings are harmless (image lacks GNU coreutils). DB persists on scratch.
- KEY RULE: the `padloc_data:/usr/local/data` overlay bind is required on EVERY padloc
  invocation (not just db-update) — padloc does the startup mkdir every time AND the DB
  lives in padloc_data. Dropping it = read-only error + no DB.

## ONE-TIME DB SETUP (must document for end users)
Both DB-backed tools need a ONE-TIME, network-connected setup run before any pipeline run
(do it on a login/data-mover node WITH internet; compute nodes have none):
  - PADLOC:        `padloc --db-update`  (with the data overlay bind) -> writes to padloc_data
  - DefenseFinder: `defense-finder update` (TODO Step 4) -> needs a writable models dir too
The wrapper/README must tell users to run this once; pipeline runs then only READ the DBs.

- 2026-06-11: DefenseFinder VERSION-MISMATCH (the dependency demon, as predicted!).
  `defense-finder update` pulled LATEST models (3.1.0 + CasFinder 3.1.1), but DF 2.0.1's
  engine only parses model format '2.0' -> `MacsypyError: ... not the right version.
  version supported is '2.0'`. Tool and DB drifted apart.
  COMPAT RULE: DF 2.0.x <-> models v2.x; models 3.x are a breaking format change needing
  newer DF. Two consistent bundles: (A) DF 2.0.1 + models 2.0.2; (B) DF 3.0.0 + models 3.1.0.
  ORIGINAL THESIS RUN (from Chp4 log) = DF **2.0.1** + models **2.0.2** -> choose Bundle A.
  Best fix: the exact 2.0.2 models already exist on Kelvin at ~/.macsyfinder/models (the
  files that made the thesis numbers) -> copy them to scratch df_models, point --models-dir
  there. No re-download, perfect fidelity. (`defense-finder update` has no version flag, so
  pinning a models version otherwise means macsydata/github-release juggling.)
- 2026-06-11: STEP 3 MILESTONE — padloc ran on pyrodigal faa/gff and called real defence
  systems (PDC-S13, target NZ_CM170101.1_60) with NO reformatting shim. The pyrodigal
  swap works for PADLOC. NOTE: padloc does NOT create its --outdir; caller must mkdir it
  first (existing run_padloc already does `mkdir -p`). Pending: compare system set vs the
  original Prokka-based padloc run on the same genome to confirm equivalence.

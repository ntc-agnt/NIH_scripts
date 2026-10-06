# smartSlurm `runAsPipeline` / `ssbatch` grammar cheat sheet

Self-contained reference for analyzing NTC bash driver scripts. Source of truth:
`/opt/ntc/smartSlurm` (`bin/runAsPipeline`, `bin/ssbatch`, `README.md`). Verified
against code, not just docs.

> **How to classify a script.** A script is a **runAsPipeline (rAP) pipeline** only
> if it contains `#@` decorator comment lines. `grep -n '#@' script.sh`. No `#@`
> lines ⇒ it is a **standard script** — read it by ordinary shell semantics, do NOT
> impose a DAG. Classification is **per script**, not per repo; one repo may hold both.

---

## 1. The `#@` job block

A *job block* = one `#@` marker line + the command line(s) directly beneath it, up to
the next **blank line**.

```
#@ stepID , dependIDs , name , reference , inputs , sbatchOptions
<command(s) to run for this step>
```

Fields are comma-split (`IFS=','`). Order and meaning:

| # | Field | Meaning | Required | Example |
|--|-------|---------|:--:|---------|
| 0 | `stepID` | Unique integer id for the step. Duplicate ⇒ run aborts. | **Yes** | `1` |
| 1 | `dependIDs` | `0` = no dependency, else upstream stepIDs **dot-joined**. | **Yes** | `0`, `1`, `1.3` |
| 2 | `name` | Program label; groups `ssbatch` resource stats + names the job. | **Yes** | `bowtie2` |
| 3 | `reference` | Ref file(s)/dir(s), **dot-joined**; rsynced to node `/tmp` under `useTmp`. | No | `genome.fa`, `db1.db2` |
| 4 | `inputs` | Input file(s), **dot-joined**; size drives resource estimation. | No | `reads.fq`, `in1.in2` |
| 5 | `sbatchOptions` | Per-step sbatch string; must start with `sbatch`. Omit ⇒ CLI default. | No | `sbatch -c 4 -t 2:0:0 --mem 24G` |

Rules verified in `runAsPipeline`:
- **Only the first three fields are required**; keep empty commas as placeholders
  (`#@2,1,merge,,,`). Empty step/dep/name ⇒ hard error.
- `sbatchOptions`, if present, **must** begin with `sbatch` or the run aborts.
- **Dot** (`.`) is the separator *inside* the dependID, reference and inputs fields
  (e.g. depends on steps 1 and 3 = `1.3`; two refs = `db1.db2`).
- `name` may carry suffixes: `name:...` (`%%:*` strips trailing `:`), `.NoCheckpoint`
  disables checkpointing for that step; otherwise the global `checkpoint` tag applies.
- **Job arrays are NOT supported** — `--array=`/`-a` in a step's sbatch options ⇒ error.
- `ssbatch` must NOT appear in a command line inside the script ⇒ error (rAP adds it).

### Where a block ends — blank line is the ONLY terminator
- A `#@` block greedily collects **every following line** (joined with `;`) until a
  **blank/whitespace-only line**. Always leave a blank line before `done`, the next
  `#@`, or any following plain command, or they get swallowed into the job command.
- Full-line `# comment`: does **not** end the block; passed through to the converted
  script but not part of the job command.
- Trailing `code # comment`: the ` #…` tail is stripped **textually** (first space-`#`
  to EOL). Avoid a literal ` #` even inside quotes (`sed 's/ #/x/'` gets truncated). A
  `#` with no preceding space (`grep '#'`) is safe.
- **Heredocs break blocks** — the parser is line-based and flattens to one line; `<<EOF`
  then swallows the rest. Don't use heredocs (or `: <<EOF` comments) inside a `#@` block.

---

## 2. The two phases (when does a line run?)

```
PHASE 1 CONVERT (once, submit host):  runAsPipeline reads script top→bottom and writes
   smartSlurmLog/slurmPipeLine.<md5>.sh   (<md5> = md5sum of your script contents)
     • for / while / #loopStart  → remember loop var (becomes part of each job flag)
     • #@ marker + command        → emit an  ssbatch --wrap "command"  call
     • any other line             → copied through unchanged
PHASE 2 SUBMIT (run the converted script):
     • plain lines  → run NOW, on the submit host
     • ssbatch call → SUBMIT a job (capture id); command runs LATER on a node once
                      -d afterok:<deps> is satisfied. Independent jobs beyond
                      firstBatchCount (5) are submitted HELD (-H), released after early
                      jobs produce ≥3 records so later jobs can be resource-estimated.
```

Consequences to watch for when reading/modifying a pipeline:
- **Plain lines that inspect a step's output run before that output exists.** Any
  `[ -f out ]`, `$(samtools view -c out)` etc. in submit scope evaluates at submission,
  before the job runs. Move such logic into a *later* `#@` step (runs on a node).
- Unchanged script ⇒ reconversion skipped (md5 match); already-successful steps skipped
  unless rerun chosen (detected by `<flag>.success` in `smartSlurmLogDir`).
- Functions defined in plain scope are **not** visible inside `#@` blocks — `export -f`
  them first.

### Two variable scopes
| Scope | What lives here | Evaluated | Visible to |
|-------|-----------------|-----------|-----------|
| **Submit** | plain lines + assignments outside any `#@` block | once, submit host (Phase 2) | later plain lines; the ref/inputs/sbatch fields + command text of every `#@` block (substituted at submit) |
| **Node** | code inside a `#@` block | later, on the compute node | only that one job |

- `$var` assigned in submit scope, used in a `#@` block ⇒ **baked in at submission**.
- `$var` assigned *and* used inside the same `#@` block ⇒ deferred, expands **on the node**.
- A var assigned inside a `#@` block exists nowhere else.
- **Steps pass data through files on shared storage, never through shell variables.**
  Convention: name each output path as a submit-scope var *immediately above* the `#@`
  that first uses it; producer writes it, consumer reads it, consumer declares the dep.
- **Submit-scope vars persist across steps, branches, loop iterations.** A missing
  assignment silently reuses a stale value — assign on *every* branch, and *inside*
  per-sample loops (not above them).

---

## 3. Loops

rAP preserves loop structure but extracts the `#@` blocks inside; the **loop variable
becomes part of each job's flag** (`1.0.findNumber.1`, `1.0.findNumber.2`, …).

- **`for` loops**: the var right after `for` is auto-detected.
  ```bash
  for file in `ls dir`; do
      #@1,0,process,,file
      process.sh $file
  done
  ```
- **`while` loops** need a hint — declare the loop var on the line above:
  ```bash
  #loopStart:f1
  while read -r f1 f2 f3 f4; do
      #@1,0,process,,f1
      process.sh "$f1"
  done < samples.txt
  ```

---

## 4. `runAsPipeline` invocation

```
runAsPipeline --script "SCRIPT [ARGS]" --tmp {useTmp|noTmp} \
   [--sbatch-options "sbatch -p short -c 1 --mem 2G -t 50:0"] \
   [--mode dryrun] [--email noEmail|noSuccEmail] \
   [--special checkpoint|excludeFailedNodes] [-A|--account=ACCT]
```

| Option | Meaning | Default |
|--------|---------|---------|
| `--script "…"` | annotated script + args, one quoted string | required |
| `--tmp useTmp\|noTmp` | copy each step's `reference` to node-local `/tmp` or not | required |
| `--sbatch-options "…"` | default sbatch for steps with no own options; must start `sbatch` | `sbatch -p short -c 1 --mem 2G -t 50:0` |
| `--mode dryrun` | build + log plan, submit nothing | submit |
| `--email noEmail\|noSuccEmail` | silence all / success emails | all |
| `--special checkpoint\|excludeFailedNodes` | checkpoint long jobs / avoid previously-failing nodes | none |

> [!WARNING]
> **Invocation-form discrepancy (verify against the DEPLOYED smartSlurm).** The
> cloned `ntc-agnt/smartSlurm` (master) `runAsPipeline` parses **named flags only**
> (`--script`, `--tmp`, …) and aborts on any unknown token. But every NTC pipeline
> README (RNAseqMapping, pipeline2, DOBBE_F) invokes the **legacy positional form**
> `runAsPipeline "SCRIPT ARGS" useTmp|noTmp run`. That positional form fails against
> this master. So the smartSlurm actually deployed at
> `/n/data1/cores/ntc/scripts/SmartSlurm` is a different version/branch, or the
> pipeline READMEs are stale. There are ≥3 smartSlurm variants in play: `ld32/SmartSlurm`
> (upstream; branch `ux-update` is the active dev line the test harness targets),
> `ntc-agnt/smartSlurm` (this clone, master), and `ntc-agnt/SmartSlurm`. **Confirm the
> deployed version and its accepted CLI before relying on either form.** TODO.

Outputs land in `smartSlurmLog/` (relative to cwd): per-step `.sh` job scripts, `.out`
logs, `.success`/`.failed` flag files, `allJobs.txt`, and `.smartSlurm.log`.
`checkRun` (run from the launch dir) browses status/logs; `cancelAllJobs` cancels this
dir's running/pending jobs.

---

## 5. `ssbatch` (smart `sbatch`) — what each step becomes

```
ssbatch [-P PROGRAM] [-I INPUTS] [-F FLAG] [SBATCH_OPTIONS] --wrap="CMD"   [dryrun]
ssbatch [-P PROGRAM] [-I INPUTS] [-F FLAG] [SBATCH_OPTIONS] SCRIPT.sh [ARGS] [dryrun]
```

- `-P` program label (stats grouping); `-I` input file(s)/dir(s) or `jobSize:N`; `-F`
  unique flag. All three optional but `-P`/`-I` make estimation good.
- **Resource estimation** (`estimateResource.sh`), needs **≥3 completed records** for a
  program/reference else uses defaults:
  - with `-I`: linear fit of input size → mem and → time from past records.
  - without `-I`: 90th percentile of the program's past mem/time.
- **jobRecord.txt**: 18-col CSV of every successful job (default `~/.smartSlurm/`).
  Estimation keys on **program (col 12)** + **reference (col 13)**, reads **memUsed
  (col 7)**, **timeUsed (col 8)**; **inputSize (col 2)** drives the size fit.
- **Partition** auto-chosen from requested run-time via `adjustPartition` in
  `config.txt` (`partition{1,2,3}TimeLimit`); a bad `-p` is corrected for you.
- **OOM/OOT auto-resubmit**: `cleanUp.sh` records usage; a job that died
  out-of-memory/time is resubmitted with **doubled** mem/time and its formula cleared.
- **Checkpointing** (opt-in): snapshot + resume a long job before it hits its limit.
- **Emails**: `cleanUp.sh` emails the job script, exact submit command, and logs (unlike
  Slurm's subject-only mail). Dampen with `noSuccEmail` / `noEmail`.
- Key `config.txt` knobs: `defaultMem=4096`M, `defaultTime=120`min, `defaultExtraMem`,
  `defaultExtraTime`, `firstBatchCount=5`.

> **`firstBatchCount` (5) ≠ the "3 records" rule.** 3 records = minimum before
> estimation works. `firstBatchCount` = how many independent rAP jobs run immediately
> before the rest are held pending early records.

---

## 6. Using `ssbatch` under other managers
`ssbatch` takes standard sbatch syntax, so Snakemake / Nextflow / Cromwell can use it as
their cluster submit command and inherit smart sizing (see smartSlurm `bin/Snakefile`,
`bin/nextflow.nf`, and README "Using ssbatch with other pipeline managers").

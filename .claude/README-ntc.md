# NTC context stub (generated — do not hand-edit)

These files are stamped into this repo by `make-repo-stub.sh` in the
`ntc-ecosystem-context` meta-repo, so a Claude Code session opened here acts as this
repo's deep-analysis agent *and* is aware of the wider NTC ecosystem:

- `.claude/skills/nih-scripts/` — this repo's own deep-analysis skill (architecture, I/O
  contract, gotchas; deeper tables under its `reference/`).
- `.claude/ntc-ecosystem.md` — the always-loaded ecosystem catalog (the other NTC tools,
  the data-flow DAG, the shared file formats). Imported by the repo's `CLAUDE.md`.
- `.claude/reference/runAsPipeline-grammar.md` — the smartSlurm grammar cheat sheet.
- `.claude/ntc-todo.md` — this repo's slice of the ecosystem BACKLOG (a generated projection;
  advisory — surface to the user, act only on confirmation).
- `.claude/ntc-changelog.md` — **repo-owned** (NOT overwritten): append a line here when you
  make a notable revision, so the ecosystem manager can track it.

Edit the originals in the meta-repo and regenerate (`make-repo-stub.sh` /
`refresh-ntc-context`); changes here will be overwritten. Generated from
ntc-ecosystem-context @ 435cddc on 2026-10-26.

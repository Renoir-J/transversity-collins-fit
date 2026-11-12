# Git and versioning policy

Keep the repository private. Use short `feat/`, `fix/`, `theory/`, `data/` or `docs/` branches and focused commits. Changes to `main` go through a pull request and squash merge. Track each deliverable in an issue with explicit acceptance evidence and a milestone.

## Local checks

```bash
git config --local core.hooksPath .githooks
git config --local pull.ff only
git config --local fetch.prune true
```

The pre-commit hook checks staged content. The pre-push hook rejects direct pushes to `main` and checks the workflow specification and repository content. Linux CI runs the complete numerical validation. Local hooks can be bypassed and do not provide server-side branch protection; verify enforcement in the repository settings.

## Versioned and local-only files

Version core Python, Julia and Fortran sources, workflow definitions, environment lock files, reviewed constraints, small synthetic fixtures and source provenance.

Keep agent instruction files, internal working notes, generated reports, run products, downloaded references, raw and normalized data caches, environments, compiled binaries and credentials local. `.gitignore` defines the exclusions. Local repository snapshots such as `config/github-state.json` are generated administrative records rather than project inputs. Example environment files may be versioned if they contain no credentials.

Adding a file to `.gitignore` does not stop tracking a previously committed file. Use `git rm --cached -- <path>` to remove it from the index while preserving the local copy. This does not remove earlier commits or change existing remote history.

## Scientific changes

Data replacements must remove the superseded likelihood. Changes to data, signs, scales, jet definitions, integration, models or evolution must invalidate affected downstream artifacts. Preserve upstream licenses and source attribution. Record run inputs, code and environment versions, seeds, checks and output hashes. Synthetic validation cannot authorize production results.

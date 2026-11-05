# Accurate Species Tree EstimatoR 2 (ASTER2❋)

**Status: alpha. ASTER2 is under active development and is *not* yet a complete replacement for ASTER.**

It is published for early testing. Command-line options, log output and result formats may change without notice between commits, and results should be cross-checked against ASTER before being used in published analyses. If you need a stable, feature-complete tool today, please use [ASTER](https://github.com/chaoszhang/ASTER).

Current pre-release: **v2.6.0** (2026-11-02).

## What is ASTER2

ASTER2 is the next-generation implementation of the ASTER package. It is written in C++20 with **no external dependency**, and a single desktop machine can handle hundreds of thousands of gene trees.

## Available now

| Executable   | Tool                         | Typical input        | Status      |
| ------------ | ---------------------------- | -------------------- | ----------- |
| `bin/astral` | ASTRAL-IV, including wASTRAL | gene trees (Newick)  | pre-release |
| `bin/caster` | CASTER-site                  | MSA (Phylip / Fasta) | pre-release |

Features in `bin/astral`:

- quartet-based species tree search with stepwise addition followed by an SPR search
- wASTRAL weighting of gene trees via `--mode 1` (hybrid, default), `2` (support), `3` (length), `4` (unweighted, equivalent to ASTRAL-III)
- branch support annotations: `sbb` (SPR block bootstrap, default), `lbb` (local block bootstrap), `qs` (quartet score)
- coalescent-unit branch lengths (`cu`, default)
- guide trees (`-g`), outgroup rooting (`--root`), multi-threading (`-t`)

Features in `bin/caster`:

- direct inference from whole-genome alignments using site patterns, without reconstructing gene trees
- handles missing taxa, multiple individuals per species, and unphased multiploid data
- extremely memory efficient (below 1 byte per SNP per sample)
- same support / branch-length annotation scheme as `bin/astral` (`sbb`, `lbb`, `qs`, `cu`)

## Not available yet

The following tools and features of ASTER are **not implemented** in ASTER2 at the moment. Please keep using [ASTER](https://github.com/chaoszhang/ASTER) for these:

- **ASTRAL-Pro / ASTRAL-Pro 3** - species tree estimation in the presence of gene duplication and loss (`astral-pro3` in ASTER).
- **CASTLES** - branch lengths in substitution-per-site units (`astral4` in ASTER integrates CASTLES-II). ASTER2 reports coalescent-unit branch lengths (`cu`) instead.
- **WASTER / WASTER-lite** - de novo phylogenomics and rapid species identification from low-coverage short reads.
- **CASTER-pair** - the site-pair flavor of CASTER; `bin/caster` here implements **CASTER-site only**.
- **dstar** - sliding-window analysis with CASTER-site.
- **Multi-locus bootstrapping** - ASTER2 does not perform bootstrap replicates. Run ASTER2 on each replicate to obtain the optimal tree, then pass these trees to ASTER with `-q`, or use the block bootstrap annotations (`sbb`, `lbb`) reported here.
- **Packaging and GUI** - no pre-built binaries (Windows / macOS), no GUI, no conda package.

## Installation

Requires `g++` version 13 or newer (C++20). No external libraries are needed.

```
make            # build
make astral     # build a single tool
make doc        # writes doc/astral.md and doc/caster.md
```

## Quick start

```
bin/astral -i example/genetree.nw -o out.stree        # from gene trees
bin/caster -i example/alignment.phylip -o out.stree   # from an alignment
```

List all options, or print the full documentation in markdown:

```
bin/astral -h    # available arguments
bin/astral -H    # full documentation
```

## Citation

If you use ASTER2, please cite the ASTER package and, when using CASTER-site, the CASTER paper:

- Chao Zhang, Rasmus Nielsen, Siavash Mirarab. ASTER: A Package for Large-scale Phylogenomic Reconstructions. *Molecular Biology and Evolution*, 2025, msaf172. <https://doi.org/10.1093/molbev/msaf172>
- Chao Zhang, Rasmus Nielsen, Siavash Mirarab. CASTER: Direct species tree inference from whole-genome alignments. *Science*, 2025. <https://doi.org/10.1126/science.adk9688>

## License

GNU Affero General Public License v3.0. See [LICENSE](LICENSE).

## Bug reports

Contact `chaozhang@pku.edu.cn`, the [ASTER user group](https://groups.google.com/forum/#!forum/aster-users), or open an issue on the [ASTER2 issues page](https://github.com/chaoszhang/aster2/issues). QQ group: 130635706.

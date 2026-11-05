# Accurate Species Tree ALgorithm (ASTRAL-IV)

ASTRAL is a tool for estimating an unrooted species tree given a set of unrooted gene trees. ASTRAL is statistically consistent under the multi-species coalescent model (and thus is useful for handling incomplete lineage sorting, i.e., ILS). ASTRAL finds the species tree that has the maximum number of shared induced quartet trees with the set of gene trees, subject to the constraint that the set of tripartitions in the species tree comes from a predefined set of tripartitions.

ASTRAL-IV re-implements [ASTRAL](https://github.com/smirarab/ASTRAL) as a scalable alternative to ASTRAL on datasets for which ASTRAL is not suitable (e.g. large datasets, multi-individual, and gene trees with missing taxa). It is written in C++20 with no external dependency, and a single desktop machine can handle hundreds of thousands of gene trees.

## Weighting modes (wASTRAL)

ASTRAL-IV generalizes the quartet score of ASTRAL-III by weighting each gene tree according to its branch support and/or branch length. Use `--mode` to choose:

* `--mode 1` (default): hybrid weighting, using both branch support and branch length
* `--mode 2`: weighting by branch support only
* `--mode 3`: weighting by branch length only
* `--mode 4`: unweighted, equivalent to the quartet score used by ASTRAL-III

Weighting down-weights gene trees that are poorly resolved or inferred from short/sparse loci, and typically improves accuracy when gene trees differ greatly in quality.

## Search

The species tree is found by stepwise addition of taxa, followed by an SPR search on the resulting tree. The number of placement rounds is controlled by `-r` (initial rounds, 8 by default) and `-s` (subsequent rounds, 8 by default):

```
bin/astral -r 16 -s 16 -i INPUT_FILE -o OUTPUT_FILE
```

For datasets with few taxa and few genes it is often worth the extra time. Guide trees can be supplied to restrict or speed up the search:

```
bin/astral -g GUIDE_TREES_IN_NEWICK_FORMAT -i INPUT_FILE -o OUTPUT_FILE
```

## Branch annotations

Internal branches can carry an annotation computed by `-u`, and branch lengths an annotation computed by `-l`:

* `sbb` - SPR block bootstrap (default for `-u`)
* `lbb` - local block bootstrap
* `qs` - quartet score
* `cu` - coalescent unit length (default for `-l`)
* `none`, or `verbose` (for `-u` only) - no annotation, or report every quantity

```
bin/astral -u lbb -l cu -i INPUT_FILE -o OUTPUT_FILE
```

## Rooting

The estimated species tree is unrooted. To root it at a given outgroup, use `--root`:

```
bin/astral --root YOUR_OUTGROUP -i INPUT_FILE -o OUTPUT_FILE
```

## Publication

Chao Zhang, Rasmus Nielsen, Siavash Mirarab, ASTER: A Package for Large-scale Phylogenomic Reconstructions, Molecular Biology and Evolution, 2025, msaf172, https://doi.org/10.1093/molbev/msaf172

## Notice

Unlike ASTRAL-III, ASTRAL-IV does not implement multi-locus bootstrapping. To obtain bootstrap supports, run ASTRAL-IV on each replicate to obtain the optimal tree for that replicate, and pass the resulting trees to ASTRAL-III using `-q`. Alternatively, use the block bootstrap annotations (`sbb` and `lbb`) reported here, which do not require re-running the search.

# Documentations
- The rest of this documentation file
- Forums (feel free to ask questions or ask for help running ASTER2):
  - [User group discussions](https://groups.google.com/forum/#!forum/aster-users)
  - [ASTER2 issues page](https://github.com/chaoszhang/aster2/issues)
  - QQ group: 130635706

## Bug Reports
Contact ``chaozhang@pku.edu.cn``, [``aster-users@googlegroups.com``](https://groups.google.com/forum/#!forum/aster-users), or post on [ASTER2 issues page](https://github.com/chaoszhang/aster2/issues).

# INSTALLATION
```
make
```

# INPUT

* The input is a single file of unrooted gene trees in Newick format, one tree per line (blank lines are ignored). See `example/genetree.nw` for an example.
* The input gene trees can have missing taxa, polytomies (unresolved branches), and multiple individuals per species. Polytomies in the input are resolved arbitrarily.
* Branch support and branch length in the input gene trees are optional; when present, they are used by the weighting modes (see the `--mode` option).
* Multiple individuals from the same species are assumed to have the same name in the input, e.g. `((name_A,name_B),(name_A,name_C));`.
* Alternatively, a mapping file can be provided using `-a`, in which case individuals mapped to the same taxon are forced to be together (monophyletic) in the species tree. The format of the mapping file is described below.

## Mapping file
ASTER2 assumes that taxon names in the input are consistent throughout.
Otherwise, it will treat them as different taxa, unless you provide a mapping file.
The mapping file should be a text file where each line contains an alias followed by a taxon name, separated by a space or tab. For example:
```
alias_A1 taxon_A
alias_A2 taxon_A
alias_B1 taxon_B
alias_B2 taxon_B
alias_B3 taxon_B
...
```
ASTER2 will internally replace all occurrences of aliases with corresponding taxon names.
If an alias in the mapping file does not appear in the input, it will be ignored.
A mapping file is also useful when you don't want to use the original taxon names in the input as taxon names in the output.
A mapping file is convient when you have multiple speciesmen per taxon and you want to specify which speciesmen belong to which taxa.

# OUTPUT
The output in is Newick format and gives:
* the species tree topology
* branch supports measured as SPR block bootstrap supports by default (>70.0 means good)
* It can also annotate branches with other quantities, such as quartet scores and bootstraps for all three topologies.

# EXECUTION
ASTER2 currently has no GUI. You need to run it through the command-line. In a terminal/PowerShell, go to the directory (location) where you have downloaded ASTER2 and issue the following command:

```
bin/astral
```

This will give you a list of options available. If you are using Windows, please replace `bin/astral` with `.\exe\astral.exe`.

To find the species tree with input from in a file called `INPUT_FILE`, use:
```
bin/astral -i INPUT_FILE
```

For example if you want to run `astral` with input `example/genetree.nw`, then run

```
bin/astral -i example/genetree.nw
```

The results will be outputted to the standard output. To save the results in a file use the `-o OUTPUT_FILE` option:
```
bin/astral -i INPUT_FILE -o OUTPUT_FILE
```

```
bin/astral -i example/genetree.nw -o example/genetree.nw.stree
```

When using a mapping file, add `-a MAPPING_FILE`:
```
bin/astral -a MAPPING_FILE -i INPUT_FILE -o OUTPUT_FILE 
```

ASTER2 supports multi-threading. To run program with 4 threads, add `-t 4`:

```
bin/astral -t 4 -i INPUT_FILE -o OUTPUT_FILE
```

ASTER2 also allows rooting at an given outgroup taxon:

```
bin/astral --root YOUR_OUTGROUP INPUT_FILE
```

# CHANGE LOG
v2.6.0 - 2026-11-02:
  * CoalescentUnitLengthCASTER: Chao Zhang - Coalescent unit length from length-weighted quartet scores (CASTER series)
  * CoalescentUnitLength: Chao Zhang - Coalescent unit branch length from quartet concordance

v2.5.2 - 2026-10-29:
  * main: Chao Zhang - Supporting annotations

v2.5.1 - 2026-10-19:
  * SPR44BlockBootstrap: Chao Zhang - memory efficiency improvements by DeepSeek
  * logQuadripartitionSupportSPR: Chao Zhang - memory efficiency improvements by DeepSeek
  * logSPR: Chao Zhang - memory efficiency improvements by DeepSeek

v2.5.0 - 2026-10-16:
  * SPR44BlockBootstrap: Yujie Wang - First version
  * logQuadripartitionSupportSPR: Yujie Wang - add quardripartition support for spr
  * logSPR: Yujie Wang - add spr

v2.4.0 - 2026-10-15:
  * Newick: Chao Zhang - AI generated Newick parser, printing error messages when throwing exceptions
  * main: Chao Zhang - Adding guide tree support

v2.3.12 - 2026-10-06:
  * BranchAnnotations: Chao Zhang - Refactoring

v2.3.11 - 2026-09-03:
  * QUADRIPARTITION_SUPPORT_POSTPROCESS: Chao Zhang - Supporting post-processing
  * BRANCH_SUPPORT: Chao Zhang - Compatible for support annotation

v2.3.10 - 2026-06-14:
  * TaxonOrderGenerator: Chao Zhang - Compatible for order prioritizing
  * optimization_algorithm::Procedure: Chao Zhang - Compatible for order prioritizing
  * TAXON_ORDER_PRIORITIZING: Chao Zhang - Compatible for order prioritizing

v2.3.9 - 2026-03-30:
  * Color: Chao Zhang - Bug fix in resolving polytomies

v2.3.8 - 2026-03-17:
  * ChangeLog: Chao Zhang - Fix negative verbose

v2.3.7 - 2026-03-02:
  * STEPWISE_COLORABLE: Chao Zhang - Signiture change of constructor and data.nElements()
  * main: Chao Zhang - Allowing no support for quad score undefined

v2.3.6 - 2026-02-26:
  * StepwiseColorPlacement: Chao Zhang - Supporting elementClearAndSetTaxonColor
  * StepwiseColorQuadripartitionScore: Chao Zhang - Supporting elementClearAndSetTaxonColor & setQuadripartitionMode
  * TwoStepPlacement: Chao Zhang - Compatible for no NNI implementation
  * optimization_algorithm::Procedure: Chao Zhang - Compatible for no NNI implementation

v2.3.5 - 2026-02-25:
  * ChangeLog: Chao Zhang - Fix formating

v2.3.4 - 2026-02-24:
  * InputParser: Chao Zhang - Add printing documentation function

v2.3.3 - 2026-02-13:
  * main: Chao Zhang - Add time in logs
  * optimization_algorithm::Procedure: Chao Zhang - Switch to sequential subsample

v2.3.2 - 2026-02-12:
  * Task: Chao Zhang - Reorder executions for better cache hits

v2.3.1 - 2026-02-08:
  * main: Chao Zhang - Add --root option

v2.3.0 - 2026-02-07:
  * LocalBlockBootstrap: Chao Zhang - First version

v2.2.0 - 2026-02-05:
  * TwoStepPlacement: Chao Zhang - Initial version

v2.1.1 - 2026-02-04:
  * StepwiseColorQuadripartitionScore: Chao Zhang - Fix interface for computing support

v2.1.0 - 2026-02-02:
  * ChangeLog: Chao Zhang - Make helper function to constructor
  * StepwiseColorNNI: Chao Zhang - Initial code
  * AnnotatedBinaryTree: Chao Zhang - Bug fix in prune above and implementation of deep copy
  * QUADRIPARTITION_STEPWISE_COLORABLE: Chao Zhang - Supporting quadripartiton
  * StepwiseColorQuadripartitionScore: Chao Zhang - Initial code
  * StepwiseColorPlacement: Chao Zhang - Refactoring for conciseness

v2.0.1 - 2026-02-01:
  * main: Chao Zhang - Change default -r -s --subsample-min
  * ChangeLog: Chao Zhang - Save date of last update
  * InputParser: Chao Zhang - Display date of latest update
  * InputParser: Chao Zhang - Exit in name conflict and ensure shortcut format
  * AnnotatedBinaryTree: Chao Zhang - Add NNI and make left heavy




# Coalescence-aware Alignment-based Species Tree EstimatoR (CASTER)

[<img src="../misc/CASTER.png" width="500"/>](../misc/CASTER.png)

Genome-wide data have the promise of dramatically improving phylogenetic inferences. Yet, inferring the true phylogeny remains a challenge, mainly because the evolutionary histories of different genomic regions differ. The traditional concatenation approach ignores such differences, resulting in both theoretical and empirical shortcomings. In response, many discordance-aware inference methods have been developed. Yet, all have their own weaknesses. Many methods rely on short recombination-free genomic segments to build gene trees and thus suffer from a lack of signals for gene tree reconstruction, resulting in poor species tree. Some methods wrongly assume that the rate of evolution is uniform across the species tree. Yet, others lack enough scalability to analyze phylogenomic data.

We introduce a new site-based species tree inference method that seeks to address these challenges without reconstructing gene trees. Our method, called CASTER (Coalescence-aware Alignment-based Species Tree EstimatoR), has two flavors: CASTER-site and CASTER-pair. The first is based on patterns in individual sites and the second is based on pairs of sites.

CASTER has several outstanding features:
1. CASTER introduces two new optimization objectives based on genomic site patterns of four species; we show that optimizing these objectives produces two estimators: CASTER-site is statistically consistent under MSC+F84 model while allowing mutation rate to change across sites and across species tree branches; CASTER-pair is statistically consistent under MSC+LM model under further assumptions.
2. CASTER comes with a scalable algorithm to optimize the objectives summed over all species quartets. Remarkably, its time complexity is linear to the number of sites and at most quasi-quadratic with respect to the number of species.
3. CASTER can handle multiple samples per species, and CASTER-site specifically can work with allele frequencies of unphased multiploid.
4. CASTER is extremenly memory efficent, requiring <1 byte per SNP per sample

Under extensive simulation of genome-wide data, including recombination, we show that both CASTER-site and CASTER-pair out-perform concatenation using RAxML-ng, as well as discordance-aware methods SVDQuartets and wASTRAL in terms of both accuracy and running time. Noticeably, CASTER-site is 60¨C150X fASTER2 than the alternative methods. It reconstructs an Avian tree of 51 species from aligned genomes with 254 million SNPs in only 3.5 hours on an 8-core desktop machine with 32 GB memory. It can also reconstruct a species tree of 201 species with approximately 2 billion SNPs using a server of 256 GB memory.

Our results suggest that CASTER-site and CASTER-pair can fulfill the need for large-scale phylogenomic inferences.

## Publication

Chao Zhang, Rasmus Nielsen, Siavash Mirarab, CASTER: Direct species tree inference from whole-genome alignments. Science (2025) https://www.science.org/doi/10.1126/science.adk9688

## Notice

Since CASTER-site and CASTER-pair assume different models, please run both and choose the result that makes more sense if you can.

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

# STOP!
Please make sure you removed paralogous alignment regions using `RepeatMasker` or alike. This will improve the accuracy of CASTER on biological datasets!

# INPUT
* The input is recommended to be a single Phylip file or vertically concatenated Phylip files in one file.
* The input can also be a single MSA in Fasta format.
* The input can also be a text file containing a list of Fasta files (one file per line).
* The input file(s) can have missing taxa and multiple speciesmen/copies per taxon.

Examples:

Single Phylip file:
```
4 3
taxon_A AAA
taxon_C CCC
taxon_G GGG
taxon_T TTT
```

Multiple Phylip files concatenated, multiploid (if using `CASTER-pair`, genes must be phased; if using `CASTER-site`, you can arbitrarily phase them):
```
6 3
taxon_A AAA
taxon_A AAA
taxon_A AAA
taxon_A AAA
taxon_C CCC
taxon_C CCC
4 2
taxon_A AA
taxon_A AA
taxon_T TT
taxon_T TT
```

Single Fasta file, multiple speciesmen:
```
>taxon_A
AAA
>taxon_A
ACA
>taxon_C
CCC
>taxon_G
GGG
>taxon_T
TTT
```
Mapping file:
```
alias_A1 taxon_A
alias_A2 taxon_A
```

Multiple Fasta file:
```
gene1.fasta
gene2.fasta
```
In `gene1.fasta`:
```
>taxon_A
AAA
>taxon_C
CCC
>taxon_G
GGG
>taxon_T
TTT
```
In `gene2.fasta` (order can change, fasta files can have missing taxa):
```
>taxon_C
CC
>taxon_A
AA
>taxon_T
TT
```

Notice: only `CASTER-site` works on unphased SNPs, you can translate VCF files into Fasta (or Phylip) in the following way.

VCF:
```
taxon_A
A/A/C/T
A/A/G/G
```
Fasta (order and phasing do not matter):
```
>taxon_A
AA
>taxon_A
AA
>taxon_A
CG
>taxon_A
TG
```

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
* branch supports measured as local block bootstrap supports by default (>95.0 means good)
* It can also annotate branches with other quantities, such as quartet scores and bootstraps for all three topologies.

# EXECUTION
ASTER2 currently has no GUI. You need to run it through the command-line. In a terminal/PowerShell, go to the directory (location) where you have downloaded ASTER2 and issue the following command:

```
bin/caster
```

This will give you a list of options available. If you are using Windows, please replace `bin/caster` with `.\exe\caster.exe`.

To find the species tree with input from in a file called `INPUT_FILE`, use:
```
bin/caster -i INPUT_FILE
```

For example if you want to run `caster` with input `example/alignment.phylip`, then run

```
bin/caster -i example/alignment.phylip
```

The results will be outputted to the standard output. To save the results in a file use the `-o OUTPUT_FILE` option:
```
bin/caster -i INPUT_FILE -o OUTPUT_FILE
```

```
bin/caster -i example/alignment.phylip -o example/alignment.phylip.stree
```

When using a mapping file, add `-a MAPPING_FILE`:
```
bin/caster -a MAPPING_FILE -i INPUT_FILE -o OUTPUT_FILE 
```

ASTER2 supports multi-threading. To run program with 4 threads, add `-t 4`:

```
bin/caster -t 4 -i INPUT_FILE -o OUTPUT_FILE
```

ASTER2 has very good parrallel efficiency up to 64 cores when input data is large. In fact, it often experiences super-linear speedup with 16 cores or more. So feel free to use as many cores as you want.

ASTER2 also allows rooting at an given outgroup taxon:

```
bin/caster --root YOUR_OUTGROUP INPUT_FILE
```

# CHANGE LOG
v2.3.5 - 2026-02-25:
  * ChangeLog: Chao Zhang - Fix formating

v2.3.4 - 2026-02-24:
  * InputParser: Chao Zhang - Add printing documentation function

v2.3.3 - 2026-02-13:
  * optimization_algorithm::Procedure: Chao Zhang - Switch to sequential subsample
  * AlignmentParser: Chao Zhang - Speeding up preprocess
  * main: Chao Zhang - Add time in logs

v2.3.2 - 2026-02-12:
  * Task: Chao Zhang - Reorder executions for better cache hits

v2.3.1 - 2026-02-08:
  * Driver: Chao Zhang - Adding more type support
  * main: Chao Zhang - Add --root option

v2.3.0 - 2026-02-07:
  * LocalBlockBootstrap: Chao Zhang - First version

v2.2.0 - 2026-02-05:
  * TwoStepPlacement: Chao Zhang - Initial version

v2.1.1 - 2026-02-04:
  * StepwiseColorQuadripartitionScore: Chao Zhang - Fix interface for computing support
  * DriverHelper: Chao Zhang - Little code refactoring, no functional change

v2.1.0 - 2026-02-02:
  * ChangeLog: Chao Zhang - Make helper function to constructor
  * AnnotatedBinaryTree: Chao Zhang - Bug fix in prune above and implementation of deep copy
  * QUADRIPARTITION_STEPWISE_COLORABLE: Chao Zhang - Supporting quadripartiton
  * StepwiseColorPlacement: Chao Zhang - Refactoring for conciseness
  * StepwiseColorQuadripartitionScore: Chao Zhang - Initial code
  * StepwiseColorNNI: Chao Zhang - Initial code
  * Color: Chao Zhang - Supporting quadripartiton

v2.0.1 - 2026-02-01:
  * main: Chao Zhang - Change default -r -s --subsample-min
  * InputParser: Chao Zhang - Display date of latest update
  * InputParser: Chao Zhang - Exit in name conflict and ensure shortcut format
  * Driver: Chao Zhang - Change prgramName to caster
  * AnnotatedBinaryTree: Chao Zhang - Add NNI and make left heavy
  * ChangeLog: Chao Zhang - Save date of last update




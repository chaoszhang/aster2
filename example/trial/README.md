# TRIAL / wTRIAL / wCASTER example

A small simulated 8-taxon alignment with a known species tree.

- `alignment.fasta` — 8 taxa, 1200 bp, no gaps
- `fasta2ref.txt` — `alignment.fasta    A` (A is the reference species)
- `true_tree.nwk` — true unrooted topology: `(((A,B),(C,D)),((E,F),(G,H)));`

Run from the repository root:

```bash
bin/trial   -i example/trial/fasta2ref.txt -o trial.tre
bin/wtrial  -i example/trial/fasta2ref.txt -o wtrial.tre
bin/wtrial  -i example/trial/fasta2ref.txt --alignment-weight -o wtrial_alignment.tre
bin/wcaster -i example/trial/fasta2ref.txt -o wcaster.tre
bin/wcaster -i example/trial/fasta2ref.txt --alignment-weight -o wcaster_alignment.tre
```

The default support annotation is `sbb`; use `-u lbb` for local block bootstrap.
All commands above should recover the true topology (RF = 0). The alignment has
strong signal (626 variable sites), so support values are high.

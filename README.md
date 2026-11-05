# Ishango Bone: supplementary analysis code

Companion code for the paper identified by DOI [10.48550/arXiv.2504.06412](https://doi.org/10.48550/arXiv.2504.06412).

## Scripts

- `ishango_analysis.py`: 50,000 sampled permutations per reference configuration; seed 42.
- `ishango_permutation_v4.py`: 1,000,000 Monte Carlo arrangements and 50,000 sampled permutations per reference configuration; seed 42.
- `ishango_reconstruction_sensitivity.py`: exhaustive enumeration of 64 permitted reconstructions; no random sampling.

The three Python files are preserved as supplied. Python 3.9 or later is recommended; only the standard library is required.

## Run

```bash
python3 ishango_analysis.py
python3 ishango_permutation_v4.py
python3 ishango_reconstruction_sensitivity.py --output-dir ishango_results
```

The permutation scripts print their results. The sensitivity script writes a text report and CSV tables, including all 64 configurations, individual and pairwise property outcomes, score distributions, choice sensitivity, matched comparisons, reference checks, and configurations satisfying all five properties.

## Property definitions

P1: G and D have equal totals.
P2: G-D row sums form a nonconstant arithmetic progression after sorting.
P3: Each M pair product equals its corresponding G-D row sum.
P4: Each of G and D admits a partition into two equal-sum pairs, regardless of their positions.
P5: Using the larger member of each M pair as its primary, the first plus fourth primary equals the second plus third.

The score is the number of properties satisfied. A complete match requires all five in the same configuration.

## Restricted reconstruction space

G is fixed at (11,13,17,19). D is either (11,21,19,9) or (21,11,19,9). The first two M pairs are either (x,6),(4,8) or (4,8),(x,6), with x in {3,4}. The third pair stays third and is (y,z), with y in {9,10} and z in {4,5}. The fourth stays fourth and is (w,7), with w in {5,4}. All six binary choices are independently combinable, giving 64 distinct configurations. The numerical edits are independent; M total is not fixed.

Exactly one configuration satisfies all five: G=(11,13,17,19), D=(21,11,19,9), M=((4,8),(4,6),(9,4),(4,7)). Individual counts are P1=64, P2=32, P3=1, P4=64, P5=16. Score counts are 2:24, 3:32, 4:7, 5:1. P1 and P4 are invariant in this restricted space.

## Interpretation and implementation notes

The 1/64 fraction is an exact descriptive proportion of the defined reconstruction space, not a p-value or historical probability. Uniqueness is conditional on the specified alternatives. P3 alone selects the complete-match configuration; the properties are not independent.

Both randomisation routines sample shuffled versions of each reference configuration's own 16-value pool. The runs labelled Monte Carlo and permutation use the same sampling mechanism. In the combined script, both random generators start at seed 42; the two streams should not be described as independent validation methods. Pools differ between reference configurations.

The original randomisation scripts display zero hits as `< 1/N`; that is a display convention, not a formal confidence bound, and zero hits do not prove zero probability. Their randomisation does not replay reconstruction edits or account automatically for exploratory property selection. Their opening comment about equal-sum pairability implying positional outer-middle equality should not be read literally: P4 is position-independent. These notes document limitations without altering the supplied scripts.

No manuscript source, figures, or new licence grant is included. Cite the paper DOI above and the specific repository commit when referencing the code.

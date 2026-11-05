#!/usr/bin/env python3
"""Exhaustive analysis of 64 permitted Ishango reconstructions.

Only the six explicitly listed binary choices vary. Counts are exact
within this configuration space, not p-values for intentional design.
P4 allows any partition into equal-sum pairs. P5 uses max(M pair).
Requires Python 3.9+; standard library only.
"""
import argparse
import csv
import json
from collections import Counter
from itertools import combinations, product
from pathlib import Path

G = (11, 13, 17, 19)
PROPERTIES = ('P1', 'P2', 'P3', 'P4', 'P5')
CHOICES = ('swap_D', 'swap_M', 'value_with_6', 'damaged_first',
           'damaged_second', 'value_with_7')
OPTIONS = ((False, True), (False, True), (3, 4), (9, 10), (4, 5), (5, 4))
LABELS = {
    'P1': 'Equal G and D totals',
    'P2': 'Nonconstant arithmetic progression of G-D row sums',
    'P3': 'M pair products equal corresponding G-D row sums',
    'P4': 'Each of G and D admits two equal-sum pairs (any pairing)',
    'P5': 'Outer and middle sums of M primaries (pair maxima) equal',
}

def is_ap(values):
    s = sorted(values)
    return s[1] > s[0] and all(
        s[i + 1] - s[i] == s[1] - s[0] for i in range(1, len(s) - 1)
    )

def has_half_pairing(column):
    return any(2 * sum(pair) == sum(column)
               for pair in combinations(column, 2))

def evaluate(g, d, m):
    row_sums = tuple(x + y for x, y in zip(g, d))
    primaries = tuple(max(pair) for pair in m)
    return {
        'P1': sum(g) == sum(d),
        'P2': is_ap(row_sums),
        'P3': all(x * y == total for (x, y), total in zip(m, row_sums)),
        'P4': has_half_pairing(g) and has_half_pairing(d),
        'P5': primaries[0] + primaries[3] == primaries[1] + primaries[2],
    }

def build_configuration(swap_D, swap_M, value_with_6, damaged_first,
                        damaged_second, value_with_7):
    d = (21, 11, 19, 9) if swap_D else (11, 21, 19, 9)
    first_two = [(value_with_6, 6), (4, 8)]
    if swap_M:
        first_two.reverse()
    m = tuple(first_two + [(damaged_first, damaged_second), (value_with_7, 7)])
    return G, d, m

def enumerate_configurations():
    rows = []
    for index, values in enumerate(product(*OPTIONS), start=1):
        g, d, m = build_configuration(*values)
        tests = evaluate(g, d, m)
        rows.append({
            'configuration_id': index,
            **dict(zip(CHOICES, values)),
            'G': g, 'D': d, 'M': m,
            'row_sums': tuple(x + y for x, y in zip(g, d)),
            'M_primaries': tuple(max(pair) for pair in m),
            'M_total': sum(sum(pair) for pair in m),
            **tests, 'score': sum(tests.values()), 'all_pass': all(tests.values()),
        })
    return rows

def reference_checks():
    references = {
        'Documented (9,4)': ((11, 21, 19, 9), ((3,6),(4,8),(9,4),(5,7))),
        'Documented (10,5)': ((11, 21, 19, 9), ((3,6),(4,8),(10,5),(5,7))),
        'Fully adjusted': ((21, 11, 19, 9), ((4,8),(4,6),(9,4),(4,7))),
    }
    checks = []
    for name, (d, m) in references.items():
        tests = evaluate(G, d, m)
        checks.append({'reference': name, **tests, 'score': sum(tests.values())})
    assert checks[-1]['score'] == 5, 'Adjusted reference must satisfy all five properties.'
    return checks

def write_csv(path, rows):
    with path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(v) if isinstance(v, (tuple, list)) else v
                             for k, v in row.items()})

def analyse(rows):
    n = len(rows)
    individual = [
        {'property': p, 'definition': LABELS[p],
         'passes': sum(r[p] for r in rows), 'configurations': n,
         'fraction': sum(r[p] for r in rows) / n}
        for p in PROPERTIES
    ]
    pairwise = [
        {'property_1': p, 'property_2': q,
         'joint_passes': sum(r[p] and r[q] for r in rows),
         'configurations': n,
         'fraction': sum(r[p] and r[q] for r in rows) / n}
        for p, q in combinations(PROPERTIES, 2)
    ]
    counts = Counter(r['score'] for r in rows)
    scores = [{'score': s, 'count': counts[s], 'fraction': counts[s] / n}
              for s in range(6)]
    conditional = []
    for choice, options in zip(CHOICES, OPTIONS):
        for value in options:
            subset = [r for r in rows if r[choice] == value]
            conditional.append({
                'choice': choice, 'value': value, 'configurations': len(subset),
                'mean_score': sum(r['score'] for r in subset) / len(subset),
                'all_pass_count': sum(r['all_pass'] for r in subset),
                **{p + '_passes': sum(r[p] for r in subset) for p in PROPERTIES},
            })
    lookup = {tuple(r[c] for c in CHOICES): r for r in rows}
    matched = []
    for i, (choice, options) in enumerate(zip(CHOICES, OPTIONS)):
        comparisons = []
        for row in rows:
            if row[choice] != options[0]:
                continue
            key = [row[c] for c in CHOICES]
            key[i] = options[1]
            comparisons.append((row, lookup[tuple(key)]))
        matched.append({
            'choice': choice, 'from_value': options[0], 'to_value': options[1],
            'matched_pairs': len(comparisons),
            'score_increases': sum(b['score'] > a['score'] for a,b in comparisons),
            'score_decreases': sum(b['score'] < a['score'] for a,b in comparisons),
            'score_unchanged': sum(b['score'] == a['score'] for a,b in comparisons),
            **{p + '_flips': sum(a[p] != b[p] for a,b in comparisons) for p in PROPERTIES},
        })
    return individual, pairwise, scores, conditional, matched

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=Path('ishango_results'))
    args = parser.parse_args(argv)
    rows = enumerate_configurations()
    assert len(rows) == 64
    assert len({(r['G'], r['D'], r['M']) for r in rows}) == 64
    checks = reference_checks()
    individual, pairwise, scores, conditional, matched = analyse(rows)
    winners = [r for r in rows if r['all_pass']]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    tables = {
        'configurations.csv': rows,
        'individual_properties.csv': individual,
        'pairwise_properties.csv': pairwise,
        'score_distribution.csv': scores,
        'choice_sensitivity.csv': conditional,
        'matched_choice_comparisons.csv': matched,
        'reference_checks.csv': checks,
    }
    for filename, table in tables.items():
        write_csv(args.output_dir / filename, table)
    if winners:
        write_csv(args.output_dir / 'all_properties_configurations.csv', winners)
    else:
        with (args.output_dir / 'all_properties_configurations.csv').open('w', newline='') as f:
            csv.writer(f).writerow(list(rows[0]))
    lines = [
        'Ishango Bone: exhaustive reconstruction sensitivity analysis',
        '64 unique configurations; six independently combinable binary choices.',
        'All configurations are counted once. No random seed or simulation is needed.',
        'Fractions describe this configuration space, not historical probabilities or p-values.',
        'P1 and P4 are invariant under these permitted changes.',
        'Changes 3->4 and 5->4 are independent; M total is not held fixed.',
        '', 'Reference checks:',
    ]
    for check in checks:
        lines.append(f"{check['reference']}: " + ', '.join(
            f'{p}={check[p]}' for p in PROPERTIES) + f"; score={check['score']}/5")
    lines.extend(['', 'Individual properties:'])
    for r in individual:
        lines.append(f"{r['property']}: {r['passes']}/64 ({r['fraction']:.6%}) - {r['definition']}")
    lines.extend(['', 'Pairwise joint counts:'])
    for r in pairwise:
        lines.append(f"{r['property_1']} & {r['property_2']}: {r['joint_passes']}/64 ({r['fraction']:.6%})")
    lines.extend(['', 'Score distribution:'])
    for r in scores:
        lines.append(f"Score {r['score']}: {r['count']}/64 ({r['fraction']:.6%})")
    lines.extend(['', f'All five properties: {len(winners)}/64 ({len(winners)/64:.6%})'])
    for r in winners:
        lines.append(f"Configuration {r['configuration_id']}: G={r['G']}; D={r['D']}; M={r['M']}")
    report = '\n'.join(lines) + '\n'
    (args.output_dir / 'report.txt').write_text(report, encoding='utf-8')
    print(report)
    print(f'CSV tables and report saved in: {args.output_dir.resolve()}')
    return rows

if __name__ == '__main__':
    main()


"""
Ishango Bone — Permutation Analysis (v4)
========================================
Changes from v3:
  - P5 (grouping rule in G and D) removed as redundant with P4
    (internal pairing): if each column subdivides into two pairs each
    summing to half the column total, the outer-inner grouping rule
    follows automatically.
  - Remaining properties renumbered to match Table 2 in the paper.

Scoring is now out of 5 properties (was 6).

Property list
-------------
P1  Equal column sums (G and D)
P2  Row sums form AP in any ordering (G and D)
P3  Column M pairs multiply to row sums
P4  Columns G and D satisfy internal pairing
P5  Grouping rule holds in Column M primaries (outer + inner pairs equal)
"""

import random
from collections import Counter

RANDOM_SEED = 42
N_PERMUTATIONS = 50_000

CONFIGS = {
    "Config 1 — original (9, 4)": {
        "G": [11, 13, 17, 19],
        "D": [11, 21, 19,  9],
        "M": [(3, 6), (4, 8), (9, 4), (5, 7)],
    },
    "Config 2 — original (10, 5)": {
        "G": [11, 13, 17, 19],
        "D": [11, 21, 19,  9],
        "M": [(3, 6), (4, 8), (10, 5), (5, 7)],
    },
    "Config 3 — fully adjusted": {
        "G": [11, 13, 17, 19],
        "D": [21, 11, 19,  9],
        "M": [(4, 8), (4, 6), (9, 4), (4, 7)],
    },
}


def is_arithmetic_progression_any_order(values):
    """
    Return True if values form a strictly arithmetic progression in any
    ordering (i.e. when sorted). Common difference must be nonzero.
    """
    s = sorted(values)
    diffs = [s[i+1] - s[i] for i in range(len(s)-1)]
    return len(set(diffs)) == 1 and diffs[0] != 0


def check_patterns(G, D, M):
    """
    Check 5 structural properties. Returns dict of name -> bool plus score.

    P1  equal_col_sums       — sum(G) == sum(D)
    P2  row_sums_AP          — four row sums form a strictly AP in any
                               ordering
    P3  mult_relation        — each Column M pair multiplies to its
                               row sum
    P4  both_cols_paired     — each column subdivides into two pairs
                               each summing to half the column total
    P5  m_grouping_rule      — outer and inner pairs of Column M
                               primaries produce equal sums
    """
    results = {}

    # P1
    sum_G = sum(G)
    sum_D = sum(D)
    results["P1_equal_col_sums"] = (sum_G == sum_D)

    # P2
    row_sums = [G[i] + D[i] for i in range(4)]
    results["P2_row_sums_AP"] = is_arithmetic_progression_any_order(row_sums)

    # P3
    results["P3_mult_relation"] = all(
        M[i][0] * M[i][1] == row_sums[i] for i in range(4)
    )

    # P4
    def has_half_pairings(col):
        s = sum(col)
        if s % 2 != 0:
            return False
        half = s // 2
        for i in range(len(col)):
            for j in range(i + 1, len(col)):
                if col[i] + col[j] == half:
                    return True
        return False

    results["P4_both_cols_paired"] = (
        has_half_pairings(G) and has_half_pairings(D)
    )

    # P5
    m_primary = [max(pair) for pair in M]
    outer_sum_m = m_primary[0] + m_primary[3]
    inner_sum_m = m_primary[1] + m_primary[2]
    results["P5_m_grouping_rule"] = (outer_sum_m == inner_sum_m)

    results["score"] = sum(
        1 for k, v in results.items() if k != "score" and v
    )
    return results


def get_pool(cfg):
    pool = list(cfg["G"]) + list(cfg["D"])
    for pair in cfg["M"]:
        pool.extend(pair)
    return pool


def rand_config(pool, rng):
    sh = pool[:]
    rng.shuffle(sh)
    G = sh[0:4]
    D = sh[4:8]
    M = [(sh[8], sh[9]), (sh[10], sh[11]),
         (sh[12], sh[13]), (sh[14], sh[15])]
    return G, D, M


def run_permutation_test(config_name, cfg, n_perms, rng):
    actual = check_patterns(cfg["G"], cfg["D"], cfg["M"])
    actual_score = actual["score"]
    pool = get_pool(cfg)

    at_least_as_good = 0
    score_dist = Counter()

    for _ in range(n_perms):
        G, D, M = rand_config(pool, rng)
        res = check_patterns(G, D, M)
        score_dist[res["score"]] += 1
        if res["score"] >= actual_score:
            at_least_as_good += 1

    return {
        "actual_results": actual,
        "at_least_as_good": at_least_as_good,
        "p_value": at_least_as_good / n_perms,
        "score_distribution": score_dist,
    }


def format_p(p):
    if p == 0:
        return f"< {1/N_PERMUTATIONS:.4%} (zero hits)"
    return f"{p:.4%}"


PATTERN_LABELS = {
    "P1_equal_col_sums":   "P1  Equal column sums (G and D)",
    "P2_row_sums_AP":      "P2  Row sums form AP in any ordering (G and D)",
    "P3_mult_relation":    "P3  Column M pairs multiply to row sums",
    "P4_both_cols_paired": "P4  Columns G and D satisfy internal pairing",
    "P5_m_grouping_rule":  "P5  Grouping rule holds in Column M primaries",
}


def print_results(config_name, cfg, test_result):
    actual = test_result["actual_results"]
    p = test_result["p_value"]
    dist = test_result["score_distribution"]
    n = N_PERMUTATIONS

    print(f"\n{'='*65}")
    print(f"  {config_name}")
    print(f"{'='*65}")
    print(f"\n  Pool  G : {cfg['G']}")
    print(f"        D : {cfg['D']}")
    print(f"        M : {cfg['M']}")

    print(f"\n  Patterns detected in actual configuration:")
    for key, label in PATTERN_LABELS.items():
        mark = "[YES]" if actual[key] else "[ no]"
        print(f"    {mark}  {label}")

    print(f"\n  Score: {actual['score']} / 5 patterns")
    print(f"\n  Permutation test ({n:,} permutations, seed={RANDOM_SEED}):")
    print(f"    Permutations scoring >= {actual['score']}: "
          f"{test_result['at_least_as_good']:,}")
    print(f"    p-value : {format_p(p)}")
    sig = "SIGNIFICANT" if p < 0.05 else "not significant"
    print(f"    Result  : {sig} at alpha = 0.05")

    print(f"\n  Score distribution across {n:,} permutations:")
    for s in sorted(dist.keys()):
        pct = dist[s] / n * 100
        bar = "#" * int(pct / 0.5)
        marker = " <-- actual" if s == actual["score"] else ""
        print(f"    Score {s}: {dist[s]:6,}  ({pct:5.2f}%)  {bar}{marker}")


def print_comparison(all_results):
    print(f"\n{'='*65}")
    print("  Cross-configuration comparison")
    print(f"{'='*65}")
    print(f"\n  {'Configuration':<35} {'Score':>7} {'p-value':>22} {'Sig?':>6}")
    print(f"  {'-'*35} {'-'*7} {'-'*22} {'-'*6}")
    for name, res in all_results.items():
        score = res["actual_results"]["score"]
        p = res["p_value"]
        sig = "yes" if p < 0.05 else "no"
        short = name.split("—")[0].strip()
        print(f"  {short:<35} {score:>4}/5  {format_p(p):>22} {sig:>6}")

    print(f"\n  Pattern-by-pattern breakdown:")
    names = list(all_results.keys())
    short_names = [
        n.split("—")[1].strip() if "—" in n else n for n in names
    ]
    header = (f"  {'Pattern':<50}"
              + "".join(f"  {s:<16}" for s in short_names))
    print(header)
    print("  " + "-" * (50 + 18 * len(names)))
    for key, label in PATTERN_LABELS.items():
        row = f"  {label:<50}"
        for res in all_results.values():
            val = "YES" if res["actual_results"][key] else "no"
            row += f"  {val:<16}"
        print(row)


def main():
    print("Ishango Bone — Permutation Analysis v4")
    print("Changes: P5 (grouping rule in G and D) removed as redundant with P4")
    print(f"Random seed : {RANDOM_SEED}")
    print(f"Permutations: {N_PERMUTATIONS:,} per configuration")

    rng = random.Random(RANDOM_SEED)
    all_results = {}

    for config_name, cfg in CONFIGS.items():
        print(f"\nRunning {config_name}...", end=" ", flush=True)
        result = run_permutation_test(
            config_name, cfg, N_PERMUTATIONS, rng
        )
        all_results[config_name] = result
        print("done.")

    for config_name, cfg in CONFIGS.items():
        print_results(config_name, cfg, all_results[config_name])

    print_comparison(all_results)


if __name__ == "__main__":
    main()


"""
Ishango Bone — Permutation Analysis (v4)
========================================
Changes from v3:
  - P5 (grouping rule in G and D) removed as redundant with P4
    (internal pairing): if each column subdivides into two pairs each
    summing to half the column total, the outer-inner grouping rule
    follows automatically.
  - Remaining properties renumbered to match Table 2 in the paper.
  - Both MC (1,000,000) and permutation test (50,000) are run.

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
N_MC = 1_000_000
N_PERM = 50_000

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
    s = sorted(values)
    diffs = [s[i+1] - s[i] for i in range(len(s)-1)]
    return len(set(diffs)) == 1 and diffs[0] != 0


def check_patterns(G, D, M):
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


def run_test(config_name, cfg, n_perms, rng):
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
        "n_perms": n_perms,
    }


def format_p(p, n):
    if p == 0:
        return f"< {1/n:.4%} (zero hits)"
    return f"{p:.4%}"


PATTERN_LABELS = {
    "P1_equal_col_sums":   "P1  Equal column sums (G and D)",
    "P2_row_sums_AP":      "P2  Row sums form AP in any ordering (G and D)",
    "P3_mult_relation":    "P3  Column M pairs multiply to row sums",
    "P4_both_cols_paired": "P4  Columns G and D satisfy internal pairing",
    "P5_m_grouping_rule":  "P5  Grouping rule holds in Column M primaries",
}


def print_results(config_name, cfg, mc_result, perm_result):
    actual = mc_result["actual_results"]
    n_mc = mc_result["n_perms"]
    n_perm = perm_result["n_perms"]
    dist = mc_result["score_distribution"]

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

    print(f"\n  Monte Carlo ({n_mc:,} arrangements, seed={RANDOM_SEED}):")
    print(f"    Hits >= {actual['score']}: {mc_result['at_least_as_good']:,}")
    print(f"    p-value : {format_p(mc_result['p_value'], n_mc)}")

    print(f"\n  Permutation test ({n_perm:,} arrangements, seed={RANDOM_SEED}):")
    print(f"    Hits >= {actual['score']}: {perm_result['at_least_as_good']:,}")
    print(f"    p-value : {format_p(perm_result['p_value'], n_perm)}")

    print(f"\n  Score distribution across {n_mc:,} MC arrangements:")
    for s in sorted(dist.keys()):
        pct = dist[s] / n_mc * 100
        bar = "#" * int(pct / 0.5)
        marker = " <-- actual" if s == actual["score"] else ""
        print(f"    Score {s}: {dist[s]:7,}  ({pct:6.3f}%)  {bar}{marker}")


def print_comparison(all_mc, all_perm):
    print(f"\n{'='*65}")
    print("  Cross-configuration comparison")
    print(f"{'='*65}")
    print(f"\n  {'Configuration':<35} {'Score':>5} {'MC hits':>10} "
          f"{'MC p-value':>22} {'Perm p-value':>22}")
    print(f"  {'-'*35} {'-'*5} {'-'*10} {'-'*22} {'-'*22}")
    for name in all_mc:
        mc = all_mc[name]
        perm = all_perm[name]
        score = mc["actual_results"]["score"]
        short = name.split("—")[0].strip()
        print(f"  {short:<35} {score:>3}/5  "
              f"{mc['at_least_as_good']:>10,}  "
              f"{format_p(mc['p_value'], mc['n_perms']):>22}  "
              f"{format_p(perm['p_value'], perm['n_perms']):>22}")

    print(f"\n  Pattern-by-pattern breakdown:")
    names = list(all_mc.keys())
    short_names = [
        n.split("—")[1].strip() if "—" in n else n for n in names
    ]
    header = (f"  {'Pattern':<50}"
              + "".join(f"  {s:<16}" for s in short_names))
    print(header)
    print("  " + "-" * (50 + 18 * len(names)))
    for key, label in PATTERN_LABELS.items():
        row = f"  {label:<50}"
        for res in all_mc.values():
            val = "YES" if res["actual_results"][key] else "no"
            row += f"  {val:<16}"
        print(row)


def main():
    print("Ishango Bone — Permutation Analysis v4")
    print("Changes: P5 (grouping rule in G and D) removed; scoring out of 5")
    print(f"Random seed      : {RANDOM_SEED}")
    print(f"MC arrangements  : {N_MC:,}")
    print(f"Perm arrangements: {N_PERM:,}")

    rng_mc   = random.Random(RANDOM_SEED)
    rng_perm = random.Random(RANDOM_SEED)

    all_mc   = {}
    all_perm = {}

    for config_name, cfg in CONFIGS.items():
        print(f"\nRunning MC   for {config_name}...", end=" ", flush=True)
        all_mc[config_name] = run_test(config_name, cfg, N_MC, rng_mc)
        print("done.")

        print(f"Running Perm for {config_name}...", end=" ", flush=True)
        all_perm[config_name] = run_test(config_name, cfg, N_PERM, rng_perm)
        print("done.")

    for config_name, cfg in CONFIGS.items():
        print_results(config_name, cfg,
                      all_mc[config_name], all_perm[config_name])

    print_comparison(all_mc, all_perm)


if __name__ == "__main__":
    main()


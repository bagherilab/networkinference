import sys
import re
import time
import pandas as pd
from itertools import groupby

# Settings.
filename = sys.argv[1]
suffix = sys.argv[2] if len(sys.argv) > 2 else ""

# Single codes.
codes_M = ["M" + str(i) for i in range(1,7)]
codes_L = ["L" + str(i) for i in range(1,7)]
codes_S = ["S" + str(i) for i in range(1,4)]
codes_N = ["N" + str(i) for i in range(1,6)]
codes_A = ["A" + str(i) for i in range(1,18)]

# Full codes.
codes_ML = [i + j for i in codes_M for j in codes_L]
codes_MLS = [i + j for i in codes_ML for j in codes_S]
codes_NS = [i + j for i in codes_N for j in codes_S]
codes_NSA = [i + j for i in codes_NS for j in codes_A]

# Job IDs.
true_ids = {x: str(y) for x, y in zip(codes_MLS, range(1, len(codes_MLS) + 1))}
null_ids = {x: str(y) for x, y in zip(codes_NSA, range(1, len(codes_NSA) + 1))}

# Algorithms.
algorithms = ["CORR", "MIDER", "GENIE3", "TIGRESS", "BANJO"]

# Dataframes.
df_missing_true = {i : pd.DataFrame(index = codes_M, columns = codes_L) for i in algorithms}
df_missing_null = {i : pd.DataFrame(index = codes_N, columns = codes_S) for i in algorithms}
df_completed_true = {i : pd.DataFrame(index = codes_M, columns = codes_L) for i in algorithms}
df_completed_null = {i : pd.DataFrame(index = codes_N, columns = codes_S) for i in algorithms}

# Regex pattens.
alg_pattern = '[ABCDEGIJMNORST3]+_'
full_pattern = 'M.L.S.'
null_pattern ='N.S.A[0-9]+'

# Process summary file for full files.
with open(filename, "r") as f:
    files = [m.group(0)[8:-4] for x in f for m in [re.search('Results_' + alg_pattern + full_pattern + suffix + '\.mat', x)] if m]
    f_algs_true = [m.group(0)[:-1] for x in files for m in [re.search(alg_pattern, x)] if m]
    f_trues = [m.group(0) for x in files for m in [re.search(full_pattern, x)] if m]

# Process summary file for null files.
with open(filename, "r") as f:
    files = [m.group(0)[8:-4] for x in f for m in [re.search('Results_' + alg_pattern + null_pattern + suffix + '\.mat', x)] if m]
    f_algs_null = [m.group(0)[:-1] for x in files for m in [re.search(alg_pattern, x)] if m]
    f_nulls = [m.group(0) for x in files for m in [re.search(null_pattern, x)] if m]

# Get missing runs for all algorithms.
for alg in algorithms:
    # Iterate through all M(otif) and L(ogic) for true models.
    alg_true_matches = [c for a, c in zip(f_algs_true, f_trues) if a == alg]
    for ML in codes_ML:
        matches = {c[4:6] for c in alg_true_matches if c[0:4] == ML}
        missing = list(set(codes_S) - matches)
        completed = ['.' if x in missing else 'x' for x in codes_S]
        df_missing_true[alg].ix[ML[0:2]][ML[2:4]] = missing
        df_completed_true[alg].ix[ML[0:2]][ML[2:4]] = ''.join(completed)

    # Iterate through all N(oise) and S(timulus) for null models.
    alg_null_matches = [c for a, c in zip(f_algs_null, f_nulls) if a == alg]
    for NS in codes_NS:
        matches = {c[4:] for c in alg_null_matches if c[0:4] == NS}
        missing = list(set(codes_A) - matches)
        completed = ['.' if x in missing else 'x' for x in codes_A]
        df_missing_null[alg].ix[NS[0:2]][NS[2:]] = missing
        df_completed_null[alg].ix[NS[0:2]][NS[2:]] = ''.join(completed)

# Print out completed report.
with open("STATUS" + suffix + ".txt", "w") as f:
    for alg in algorithms:
        f.write('== ' + alg + ' =' + '='*(27 - len(alg)) + '\n\n')
        f.write(df_completed_true[alg].to_string() + '\n\n')
        f.write(df_completed_null[alg].to_string() + '\n\n')

# Calculates continuous job ids.
def get_job_ids(ids):
    sort = sorted([int(x) for x in ids])
    ranges = []
    for key, group in groupby(enumerate(sort), lambda x: x[0] - x[1]):
        g = list(group);
        ranges.append((g[0][1], g[-1][1]))
    return ''.join(['[' + str(i) + '-' + str(j) + '] (' + str(j - i + 1) + ')\n' for i, j in ranges])

# Print out missing report.
with open("MISSING" + suffix + ".txt", "w") as f:
    for alg in algorithms:
        f.write('== ' + alg + ' =' + '='*(20 - len(alg)) + '\n\n')
        true_jobs = []
        null_jobs = []

        for ML in codes_ML:
            missing = df_missing_true[alg].ix[ML[0:2]][ML[2:4]]
            for miss in missing:
                code = ML + miss
                true_jobs.append(true_ids[code])
        f.write(get_job_ids(true_jobs) + '\n')

        f.write("-" * 25 + "\n\n")

        for NS in codes_NS:
            missing = df_missing_null[alg].ix[NS[0:2]][NS[2:]]
            for miss in missing:
                code = NS + miss
                null_jobs.append(null_ids[code])
        f.write(get_job_ids(null_jobs) + '\n')

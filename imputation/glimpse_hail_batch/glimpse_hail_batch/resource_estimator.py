from argparse import ArgumentParser
from collections import namedtuple
from tabulate import tabulate
import statistics
import math
import re

from .globals import find_chunks


ChunkStats = namedtuple('ChunkStats', ['batch_size', 'mean_est_cpu', 'max_est_cpu', 'min_est_cpu', 'mean_est_memory', 'max_est_memory', 'min_est_memory'])


def estimate_resources(args: dict):
    chunks = find_chunks(args['reference_dir'],
                         args['chunk_info_dir'],
                         re.compile(args['binary_reference_file_regex']),
                         requester_pays_config=args['gcs_requester_pays_configuration'])

    chunks = [chunk for chunk in chunks]
    chunk_sizes = [chunk.n_rare + chunk.n_common for chunk in chunks]

    mean_chunk_size = statistics.mean(chunk_sizes)
    max_chunk_size = max(chunk_sizes)

    goal_runtime = args['target_runtime_mins']
    max_runtime = args['max_runtime_mins']

    results_table = []

    for cores in (1, 2, 4, 8, 16):
        max_batch_size =  cores * max_runtime / (1e-5 * max_chunk_size)  # 5e-6
        optimal_batch_size = cores * goal_runtime / (1e-5 * mean_chunk_size)
        actual_batch_size = min(optimal_batch_size, max_batch_size)
        n_batches = math.ceil(args['n_samples']  / actual_batch_size)

        max_memory = 0
        mean_memory = 0
        for chunk in chunks:
            max_memory = max(max_memory, math.ceil((800e-3 + 0.97e-6 * chunk.n_rare * cores + 14.6e-6 * chunk.n_common * cores + 6.5e-9 * (chunk.n_rare + chunk.n_common) * max_batch_size + 13.7e-3 * max_batch_size + 1.8e-6 * (chunk.n_rare + chunk.n_common) * math.log(max_batch_size)))) / cores
            mean_memory = max(mean_memory, math.ceil((800e-3 + 0.97e-6 * chunk.n_rare * cores + 14.6e-6 * chunk.n_common * cores + 6.5e-9 * (chunk.n_rare + chunk.n_common) * optimal_batch_size + 13.7e-3 * optimal_batch_size + 1.8e-6 * (chunk.n_rare + chunk.n_common) * math.log(optimal_batch_size)))) / cores

        rough_cost = n_batches * cores * (goal_runtime / 60) * 0.03 * len(chunks)
        rough_cost_per_sample = rough_cost / args['n_samples']

        results_table.append([cores, args['n_samples'], goal_runtime, max_runtime, optimal_batch_size, max_batch_size, actual_batch_size, max_memory, mean_memory, rough_cost, rough_cost_per_sample])

    print(tabulate(results_table, headers=['cores', 'n_samples', 'goal_runtime', 'max_runtime', 'optimal_batch_size', 'max_batch_size', 'actual_batch_size', 'max_memory', 'mean_memory', 'rough_cost_estimate', 'rough_cost_estimate_per_sample']))


if __name__ == '__main__':
    p = ArgumentParser()

    p.add_argument("--reference-dir", type=str, required=True)
    p.add_argument("--chunk-info-dir", type=str, required=True)
    p.add_argument('--binary-reference-file-regex', type=str, required=True)

    p.add_argument('--max-runtime-mins', type=int, required=True)
    p.add_argument('--target-runtime-mins', type=int, required=True)
    p.add_argument('--n-samples', type=int, required=True)

    p.add_argument('--gcs-requester-pays-configuration', type=str, required=False)

    args = p.parse_args()

    estimate_resources(vars(args))

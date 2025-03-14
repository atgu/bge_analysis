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
                         re.compile(args['chunk_file_regex']),
                         requester_pays_config=args['gcs_requester_pays_configuration'])

    results_table = []

    for actual_batch_size in range(50, 2000, 50):
        for cores in (1, 2, 4, 8, 16):
            n_batches = math.ceil(args['n_samples']  / actual_batch_size)

            rough_cost = 0

            n_lowmem = 0
            n_standard = 0
            n_highmem = 0
            n_oom = 0

            max_runtime = 0
            min_runtime = 1000000000
            runtimes = []

            for chunk in chunks:
                memory_req = math.ceil((800e-3 + 0.97e-6 * chunk.n_rare * cores + 14.6e-6 * chunk.n_common * cores + 6.5e-9 * (chunk.n_rare + chunk.n_common) * actual_batch_size + 13.7e-3 * actual_batch_size + 1.8e-6 * (chunk.n_rare + chunk.n_common) * math.log(actual_batch_size))) / cores

                if memory_req <= 1:
                    cost_per_core_hour = 0.02429905
                    n_lowmem += 1
                elif 1 < memory_req <= 4:
                    cost_per_core_hour = 0.02684125
                    n_standard += 1
                elif 4 <= memory_req < 7.5:
                    cost_per_core_hour = 0.02929425
                    n_highmem += 1
                else:
                    cost_per_core_hour = 0.02929425
                    n_oom += 1

                estimated_runtime_mins = 3 + ((1e-5 * (chunk.n_rare + chunk.n_common) * actual_batch_size) / cores)

                rough_cost += n_batches * cores * (estimated_runtime_mins / 60) * cost_per_core_hour

                max_runtime = max(max_runtime, estimated_runtime_mins)
                min_runtime = min(min_runtime, estimated_runtime_mins)
                runtimes.append(estimated_runtime_mins)

            rough_cost_per_sample = rough_cost / args['n_samples']
            mean_runtime = statistics.mean(runtimes)

            results_table.append([cores, args['n_samples'], n_batches, min_runtime, max_runtime, mean_runtime, actual_batch_size, rough_cost, rough_cost_per_sample, n_lowmem, n_standard, n_highmem, n_oom])

    print(tabulate(results_table, headers=['cores', 'n_samples', 'n_batches', 'min_runtime', 'max_runtime', 'mean_runtime', 'batch_size', 'rough_cost_estimate', 'rough_cost_estimate_per_sample', 'n_lowmem', 'n_standard', 'n_highmem', 'n_oom']))


if __name__ == '__main__':
    p = ArgumentParser()

    p.add_argument("--reference-dir", type=str, required=True)
    p.add_argument("--chunk-info-dir", type=str, required=True)
    p.add_argument('--binary-reference-file-regex', type=str, required=True)
    p.add_argument('--chunk-file-regex', type=str, required=True)
    p.add_argument('--n-samples', type=int, required=True)
    p.add_argument('--gcs-requester-pays-configuration', type=str, required=False)

    args = p.parse_args()

    estimate_resources(vars(args))

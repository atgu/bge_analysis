import pandas as pd
import io
import re
import json
import os
import uuid
from collections import defaultdict, namedtuple
from humanize import naturalsize
from typing import Dict, List, Optional, Tuple, Union

import hailtop.fs as hfs


Chunk = namedtuple('Chunk', ['contig', 'chunk_idx', 'n_common', 'n_rare', 'info_path', 'path', 'n_variants'])


def reference_file_dir_str(output_dir: str) -> str:
    output_dir = output_dir.rstrip('/')
    return f'{output_dir}/binary_reference/'


def chunk_info_file_dir_str(output_dir: str) -> str:
    output_dir = output_dir.rstrip('/')
    return f'{output_dir}/chunks/'


def chunks_file_str(chunks_dir: str, contig: str) -> str:
    chunks_dir = chunks_dir.rstrip('/')
    return f'{chunks_dir}/chunks_{contig}.txt'


def reference_file_str() -> str:
    return f'ref_chunk_${{CONTIG}}_${{CHUNKINDEX}}'


def parse_reference_chunk_path(path: str) -> Tuple[str, int]:
    segments = path.split('_')
    contig = '_'.join(segments[2:-1])
    chunk_index = int(segments[-1])
    return (contig, chunk_index)


chunk_manifest_header = ['chunk_idx', 'contig', 'buffered_region', 'actual_region', 'cm_dummy', 'mb_dummy', 'n_variants', 'n_common']


def find_chunks(split_reference_output_dir: str,
                *,
                requested_contig: Optional[str] = None,
                requested_chunk_index: Optional[int] = None) -> List[Chunk]:
    chunks = []
    chunk_info_dir = chunk_info_file_dir_str(split_reference_output_dir)
    reference_dir = reference_file_dir_str(split_reference_output_dir)

    reference_files = hfs.ls(reference_dir)
    chunk_info_files = hfs.ls(chunk_info_dir)

    reference_paths = {}
    for file in reference_files:
        contig, chunk_index = parse_reference_chunk_path(file.path[len(reference_dir):])
        reference_paths[(contig, chunk_index)] = file.path

    for chunk_info_file in chunk_info_files:
        with hfs.open(chunk_info_file.path, 'r') as f:
            chunks_df = pd.read_csv(io.StringIO(f.read()), sep="\t", names=chunk_manifest_header)
        for _, row in chunks_df.iterrows():
            chunks.append(Chunk(contig=row['contig'],
                                chunk_idx=int(row['chunk_idx']),
                                n_common=row['n_common'],
                                n_rare=row['n_variants'] - row['n_common'],
                                info_path=chunk_info_file.path,
                                path=reference_paths[(row['contig'], int(row['chunk_idx']))],
                                n_variants=row['n_variants']))

    chunks.sort(key=lambda c: c.chunk_idx)

    if requested_contig is not None:
        chunks = [chunk for chunk in chunks if chunk.contig == requested_contig]
        if requested_chunk_index is not None:
            chunks = [chunk for chunk in chunks if chunk.chunk_idx == requested_chunk_index]

    return chunks


def get_bucket(path: str):
    return path.split('/')[2]


def rewrite_path(mount_path: str, path: str):
    assert mount_path.endswith('/')
    file_path = '/'.join(path.split('/')[3:])
    return f'{mount_path}{file_path}'


class Sample:
    @staticmethod
    def from_json(d: dict) -> 'Sample':
        return Sample(d['sample_id'], d['cram_path'], d['cram_index_path'])

    def __init__(self, sample_id: str, cram_path: str, cram_index_path: str):
        self.sample_id = sample_id
        self.cram_path = cram_path
        self.cram_index_path = cram_index_path
        self.cram_basename = os.path.basename(self.cram_path)
        self.cram_index_basename = os.path.basename(self.cram_index_path)
        self.uuid = uuid.uuid4().hex[:8]

    def remote_cram_path(self, cram_temp_dir: str) -> str:
        cram_temp_dir = cram_temp_dir.rstrip('/')
        return f'{cram_temp_dir}/{self.cram_basename}'

    def remote_cram_index_path(self, cram_temp_dir: str) -> str:
        cram_temp_dir = cram_temp_dir.rstrip('/')
        return f'{cram_temp_dir}/{self.cram_index_basename}'

    def local_cram_path(self, mount_path: str, cram_temp_dir: str) -> str:
        remote_cram_path = self.remote_cram_path(cram_temp_dir)
        return rewrite_path(mount_path, remote_cram_path)

    def local_cram_index_path(self, mount_path: str, cram_temp_dir: str) -> str:
        remote_crai_path = self.remote_cram_index_path(cram_temp_dir)
        return rewrite_path(mount_path, remote_crai_path)

    def to_dict(self):
        return {'sample_id': self.sample_id, 'cram_path': self.cram_path, 'cram_index_path': self.cram_index_path}


class SampleGroup:
    @staticmethod
    def from_json(d: dict) -> 'SampleGroup':
        samples = [Sample.from_json(s) for s in d['samples']]
        return SampleGroup(samples, d['sample_group_index'], d['output_dir'], d['mount_path'])

    def __init__(self,
                 samples: List[Sample],
                 sample_group_index: int,
                 output_dir: str,
                 mount_path: str):
        self.samples = samples
        self.sample_group_index = sample_group_index
        self._output_dir = output_dir.rstrip('/')
        self.mount_path = mount_path
        self.uuid = uuid.uuid4().hex[:8]

        sample_base_names = [os.path.basename(sample.cram_path) for sample in samples]
        assert len(set(sample_base_names)) == len(sample_base_names)

    def to_dict(self):
        samples = [s.to_dict() for s in self.samples]
        return {
            'samples': samples,
            'sample_group_index': self.sample_group_index,
            'output_dir': self._output_dir,
            'mount_path': self.mount_path,
        }

    @property
    def temp_dir(self):
        return f'{self._output_dir}/sg-{self.sample_group_index}'

    @property
    def remote_cram_temp_dir(self):
        return f'{self.temp_dir}/crams/'

    @property
    def cram_list_output_file(self):
        return f'{self.temp_dir}/phase/crams.list'

    def cram_list_gcloud_output_file(self, start: int, end: int):
        return f'{self.temp_dir}/copy-crams/{start}-{end}/gcloud_crams.list'

    def phased_output_file_root(self, contig: str, chunk_index: int):
        return f'{self.temp_dir}/phase/{contig}/chunk-{chunk_index}'

    def ligate_output_file_root(self, contig: str):
        return f'{self.ligate_output_dir}{contig}'

    @property
    def ligate_output_dir(self):
        return f'{self.temp_dir}/ligate/'

    @property
    def merged_mt_path(self):
        return f'{self.temp_dir}/merge-vcf/merged.mt'

    @property
    def success_file(self):
        return f'{self.temp_dir}/SUCCESS'

    @property
    def coverage_metric_file(self):
        return f'{self.temp_dir}/coverage_metrics.txt.gz'

    def phased_glimpse_checkpoint_file(self, contig: str, chunk_index: int):
        return f'{self.temp_dir}/phase/{contig}/chunk-{chunk_index}/glimpse_checkpoint'

    def initialize_phased_glimpse_checkpoint_file(self, contig: str, chunk_index: int) -> str:
        path = self.phased_glimpse_checkpoint_file(contig, chunk_index)
        with hfs.open(path, 'w') as f:
            f.write("")
        return path

    def get_phased_output_file_names(self, contig_chunks: Dict[str, List[Chunk]]) -> Dict[str, List[str]]:
        output_files = defaultdict(list)
        for contig, chunks in contig_chunks.items():
            for chunk in chunks:
                output_files[contig].append(self.phased_output_file_root(contig, chunk.chunk_idx))
        return output_files

    def get_ligate_output_file_names(self, contig_chunks: Dict[str, List[Chunk]]) -> Dict[str, str]:
        output_files_by_contig = {}
        for contig in contig_chunks.keys():
            output_files_by_contig[contig] = self.ligate_output_file_root(contig)
        return output_files_by_contig

    @property
    def local_crams(self):
        return [sample.local_cram_path(self.mount_path, self.remote_cram_temp_dir) for sample in self.samples]

    @property
    def local_cram_indices(self):
        return [sample.local_cram_index_path(self.mount_path, self.remote_cram_temp_dir) for sample in self.samples]

    @property
    def remote_crams(self):
        return [sample.remote_cram_path(self.remote_cram_temp_dir) for sample in self.samples]

    @property
    def remote_cram_indices(self):
        return [sample.remote_cram_index_path(self.remote_cram_temp_dir) for sample in self.samples]

    @property
    def original_crams(self):
        return [sample.cram_path for sample in self.samples]

    @property
    def original_cram_indices(self):
        return [sample.cram_index_path for sample in self.samples]

    @property
    def sample_ids(self):
        return [sample.sample_id for sample in self.samples]

    def write_cram_list(self) -> str:
        with hfs.open(self.cram_list_output_file, 'w') as f:
            for cram, cram_idx, sample_id in zip(self.local_crams, self.local_cram_indices, self.sample_ids):
                f.write(f'{cram}##idx##{cram_idx} {sample_id}\n')
        return self.cram_list_output_file

    def write_gcloud_cram_copy_list(self, start: int, end: int) -> str:
        src_files = self.original_crams[start:end] + self.original_cram_indices[start:end]
        output_file = self.cram_list_gcloud_output_file(start, end)
        with hfs.open(output_file, 'w') as f:
            for file in src_files:
                f.write(file + '\n')
        return output_file

    def write_sample_group_dict(self):
        with hfs.open(self.temp_dir + '/sample_group.json', 'w') as f:
            f.write(json.dumps(self.to_dict(), indent=4) + '\n')

    @property
    def name(self):
        return f'sample-group-{self.sample_group_index}'


def find_crams(sample_manifest: str, sample_id_col: str, cram_path_col: str, cram_index_path_col: str, n_samples: Optional[int]) -> List[Sample]:
    with hfs.open(sample_manifest, 'r') as f:
        manifest = pd.read_csv(io.StringIO(f.read()), sep='\t')
    sample_ids = manifest[sample_id_col].to_list()
    cram_paths = manifest[cram_path_col].to_list()
    cram_index_paths = manifest[cram_index_path_col].to_list()

    samples = [Sample(sample_id, cram, cram_index) for sample_id, cram, cram_index in zip(sample_ids, cram_paths, cram_index_paths)]
    if n_samples is not None:
        samples = samples[:n_samples]

    return samples


def split_samples_into_groups(samples: List[Sample],
                              desired_group_size: int,
                              output_dir: str,
                              mount_path: str) -> List[SampleGroup]:
    groups: List[SampleGroup] = []
    working_group = []
    sample_group_idx = 0
    for sample in samples:
        if len(working_group) >= desired_group_size:
            groups.append(SampleGroup(working_group, sample_group_idx, output_dir, mount_path))
            working_group = []
            sample_group_idx += 1
        working_group.append(sample)

    groups.append(SampleGroup(working_group, sample_group_idx, output_dir, mount_path))

    return groups


def get_ligate_storage_requirement(base_gib: int, n_samples: int, n_variants: int) -> str:
    n_bytes = (base_gib * 1024 * 1024 * 1024) + int(0.91527225871 * n_samples * n_variants)
    size = naturalsize(n_bytes, binary=True, format="%d")
    return size.replace(' ', '')

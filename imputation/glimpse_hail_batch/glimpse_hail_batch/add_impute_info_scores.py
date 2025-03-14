import argparse
import pandas as pd
from typing import List

import hail as hl
from hail.matrixtable import MatrixTable


def N(mt: MatrixTable):
    return hl.agg.count_where(hl.is_defined(mt.GP))


def e_i_j(mt: MatrixTable):
    return mt.GP[1] + 2 * mt.GP[2]


def f_i_j(mt: MatrixTable):
    return mt.GP[1] + 4 * mt.GP[2]


def theta_hat(mt: MatrixTable):
    return hl.agg.sum(e_i_j(mt)) / (2 * N(mt))


def AF(mt: MatrixTable):
    return hl.agg.sum(mt.GT.n_alt_alleles()) / (2 * N(mt))


def IMPUTE_INFO(mt: MatrixTable):
    return hl.if_else((theta_hat(mt) == 0) | (theta_hat(mt) == 1),
                      1,
                      1 - hl.agg.sum(f_i_j(mt) - e_i_j(mt) ** 2) / (2 * N(mt) * theta_hat(mt) * (1 - theta_hat(mt))))


def main(input_path: str,
         output_path: str,
         sample_annotations: str,
         sample_id_col_input: str,
         sample_id_col_ann: str,
         group_by_cols: List[str],
         overwrite: bool,
         init_kwargs: dict):
    assert output_path.endswith('.mt')

    hl.init(**init_kwargs)

    pd_df = pd.read_csv(sample_annotations, sep="\t", usecols=[sample_id_col_ann, *group_by_cols], dtype=str)
    sample_ann = hl.Table.from_pandas(pd_df, sample_id_col_ann)

    if input_path.endswith('.vcf.bgz'):
        mt = hl.import_vcf(input_path, reference_genome="GRCh38")
    else:
        mt = hl.read_matrix_table(input_path)

    mt = mt.annotate_cols(**sample_ann[mt[sample_id_col_input]])

    mt = mt.annotate_rows(info=mt.info.annotate(impute=hl.struct(AF=AF(mt), INFO=IMPUTE_INFO(mt))))

    for col in group_by_cols:
        mt = mt.annotate_rows(info=mt.info.annotate(**{col: hl.agg.group_by(mt[col], hl.struct(INFO=IMPUTE_INFO(mt), AF=AF(mt)))}))

    mt.write(output_path, overwrite=overwrite)
    mt_count = hl.read_matrix_table(output_path)

    print(mt_count.count())
    print(mt.describe())


if __name__ == '__main__':
    class ParseKeyValuePairs(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            extra_args = {}
            for value in values:
                try:
                    key, val = value.split("=", 1)
                    extra_args[key] = val
                except ValueError:
                    parser.error(f"Invalid format '{value}'. Expected key=value.")
            setattr(namespace, self.dest, extra_args)

    parser = argparse.ArgumentParser()

    parser.add_argument('--input-path', required=True, type=str)
    parser.add_argument('--output-path', required=True, type=str)
    parser.add_argument('--sample-annotations', required=True, type=str)
    parser.add_argument('--sample-id-col-input', required=False, type=str, default='s')
    parser.add_argument('--sample-id-col-ann', required=True, type=str)
    parser.add_argument('--group-by-col', type=str, required=True, nargs='+')
    parser.add_argument('--overwrite', action='store_true', required=False)
    parser.add_argument('--hl-init-kwarg', action=ParseKeyValuePairs, nargs='*', dest='init_kwargs')

    args = parser.parse_args()

    main(args.input_path,
         args.output_path,
         args.sample_annotations,
         args.sample_id_col_input,
         args.sample_id_col_ann,
         args.group_by_col,
         args.overwrite,
         args.init_kwargs)

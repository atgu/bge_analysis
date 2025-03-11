#!/usr/bin/env python3

from setuptools import setup, find_packages


setup(
    name='glimpse_hail_batch',
    author="Jackie Goldstein",
    author_email="jigold@broadinstitute.org",
    description="GLIMPSE pipeline implemented in Hail Batch.",
    packages=find_packages(),
    package_dir={'glimpse_hail_batch': 'glimpse_hail_batch'},
    package_data={"glimpse_hail_batch": ["data/*"]},
    python_requires=">=3.9",
    include_package_data=True,
    install_requires=[
        'hail',
        'pandas>=2,<3',
    ],
)

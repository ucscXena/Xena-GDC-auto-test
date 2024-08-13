import argparse
import logging
import time
import os
import subprocess
import sys

import geneExpressionValidation as ge_test
import mirnaValidation as mirna_test
import cnvSegmentedValidation as seg_cnv_test
import cnvGeneLevelValidation as gl_cnv_test
import somaticMutationValidation as somaticmutation_test
import methylationValidation as methylation_test


valid_dtype = [
    'star_counts',
    'star_tpm',
    'star_fpkm',
    'star_fpkm-uq',
    'mirna',
    'mirna_isoform',
    'segment_cnv_ascat-ngs', 
    'segment_cnv_DNAcopy',
    'masked_cnv_DNAcopy',
    'allele_cnv_ascat2',
    'allele_cnv_ascat3',
    'somaticmutation_wxs',
    'somaticmutation_targeted',
    'gene-level_ascat-ngs',
    'gene-level_ascat2',
    'gene-level_ascat3',
    'gene-level_absolute',
    'methylation_epic',
    'methylation_epic_v2',
    'methylation27',
    'methylation450',
    'protein',
    'clinical',
    'survival',
]


def create_parser():
    """
    Construct the program options.
    """

    parser = argparse.ArgumentParser(
        description='Run tests on Xena matrices',
    )
    parser.add_argument(
        '-p',
        '--projects',
        type=str,
        nargs='+',
        required=True,
        help='The project name.',
    )
    parser.add_argument(
        '-t',
        '--datatype',
        type=str,
        nargs='+',
        required=False,
        help='The Xena data type of the file.',
    )

    return parser


def run_tests(project, data_type):
    """Run tests on Xena matrices.

    Args:
        project (str): The project to have tests run on. 

    Returns:
        result (str): 'PASSED' or 'FAILED' for test.
    """

    path = '../{}/Xena_Matrices/{}.{}.tsv'.format(project, project, data_type)
    if data_type.startswith('star_'):
        result = ge_test.main(data_type, path, project)
    elif data_type == 'mirna':
        result = mirna_test.main(project, path, data_type)
    elif data_type == 'segment_cnv_ascat-ngs' or data_type == 'masked_cnv_DNAcopy':
        result = seg_cnv_test.main(project, path, data_type)
    elif data_type == 'gene-level_ascat-ngs' or data_type == 'gene-level_ascat2' or data_type == 'gene-level_ascat3' or data_type == 'gene-level_absolute':
        result = gl_cnv_test.main(project, path, data_type)
    elif data_type == 'somaticmutation_wxs' or data_type == 'somaticmutation_targeted': 
        result = somaticmutation_test.main(project, path, data_type)
    elif data_type == 'methylation_epic' or data_type == 'methylation_epic_v2' or data_type == 'methylation27' or data_type == 'methylation450':
        result = methylation_test.main(project, path, data_type)
    # elif data_type == 'protein':
    # elif data_type == 'survival':
    # elif data_type == 'allele_cnv_ascat2':
    # elif data_type == 'allele_cnv_ascat3':
        
    return result


def main():
    logger = logging.getLogger(__name__)
    handlers = logging.StreamHandler(sys.stdout), logging.FileHandler(os.path.join('test_' + time.strftime("%Y%m%d-%H%M%S") + '.log',))
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        handlers=handlers,
    )
    parser = create_parser()
    options = parser.parse_args()
    test_results = []
    projects = options.projects
    data_types = options.datatype
    if data_types is not None:
        for dt in data_types:
             if dt not in valid_dtype:
                logger.info('ValueError: Unsupported data type: {}.'.format(dt))
                exit(1)
    for project in projects:
        if data_types is None:
            data_types = os.listdir('../{}/Raw_Data/'.format(project))
        for data_type in data_types:
            if data_type == 'STAR':
                for star_dtype in ['star_counts', 'star_tpm', 'star_fpkm', 'star_fpkm-uq']:
                    result = run_tests(project, star_dtype)
                    test_results.append([project, star_dtype, result])
            elif data_type == 'clinical' or data_type == 'survival' or data_type == 'protein' or data_type =='allele_cnv_ascat2' or data_type == 'allele_cnv_ascat3':
                logger.info('{} test still in progress of being written or having bugs fixed.'.format(data_type))
            else: 
                result = run_tests(project, data_type)
                test_results.append([project, data_type, result])
        data_types = None
    for r in test_results:
        logger.info('{} data for {} has {}.'.format(r[1], r[0], r[2]))

if __name__ == '__main__':
    main()

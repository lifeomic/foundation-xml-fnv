from unittest import TestCase
from src.convert import extract_fusion_variant
from src.convert import gather_attributes
from src.convert import calculate_interpretation

expected_results_no_other_gene = {
    'FusionVariants': [{
        'sample_id': 'SA-1612348',
        'gene1': 'NF1',
        'gene2': 'N/A',
        'effect': 'truncation',
        'chromosome1': 'chr17',
        'start_position1': 29557687,
        'end_position1': 29887856,
        'chromosome2': 'chr6',
        'start_position2': 66426718,
        'end_position2': 66427149,
        'interpretation': 'Likely pathogenic',
        'sequence_type': 'somatic',
        'in-frame': 'unknown',
        'attributes': {
            'equivocal': 'true',
            'supporting-read-pairs': 83
        }
    }]
}

foundation_source_dict_no_other_gene = {
    'variant-report': {
        'rearrangements': {
            'rearrangement': [
                {
                    '@equivocal': 'true',
                    '@in-frame': 'unknown',
                    '@other-gene': 'N/A',
                    '@pos1': 'ch17:29557687-29887856',
                    '@pos2': 'ch6:66426718-66427149',
                    '@status': 'likely',
                    '@supporting-read-pairs': 83,
                    '@targeted-gene': 'NF1',
                    '@type': 'truncation',
                    'dna-evidence': {
                        '@sample': 'SA-1612348'
                    }
                }
            ]
        }
    }
}

expected_results_with_other_gene = {
    'FusionVariants': [{
        'sample_id': 'SA-1612348',
        'gene1': 'NF1',
        'gene2': 'PIM1',
        'effect': 'truncation',
        'chromosome1': 'chr17',
        'start_position1': 29557687,
        'end_position1': 29887856,
        'chromosome2': 'chr6',
        'start_position2': 66426718,
        'end_position2': 66427149,
        'interpretation': 'Likely pathogenic',
        'sequence_type': 'somatic',
        'in-frame': 'unknown',
        'attributes': {
            'equivocal': 'true',
            'supporting-read-pairs': 83
        }
    }]
}


foundation_source_dict_with_other_gene = {
    'variant-report': {
        'rearrangements': {
            'rearrangement': [
                {
                    '@equivocal': 'true',
                    '@in-frame': 'unknown',
                    '@other-gene': 'PIM1',
                    '@pos1': 'ch17:29557687-29887856',
                    '@pos2': 'ch6:66426718-66427149',
                    '@status': 'likely',
                    '@supporting-read-pairs': 83,
                    '@targeted-gene': 'NF1',
                    '@type': 'truncation',
                    'dna-evidence': {
                        '@sample': 'SA-1612348'
                    }
                }
            ]
        }
    }
}


class ConvertTest(TestCase):

    def test_convert_no_other_gene(self):
        fnv_resources = extract_fusion_variant(foundation_source_dict_no_other_gene)
        print(fnv_resources)
        self.maxDiff = None
        self.assertDictEqual(expected_results_no_other_gene, fnv_resources)

    def test_convert_with_other_gene(self):
        fnv_resources = extract_fusion_variant(foundation_source_dict_with_other_gene)
        print(fnv_resources)
        self.maxDiff = None
        self.assertDictEqual(expected_results_with_other_gene, fnv_resources)

    def test_calculate_interpretation(self):
        self.assertEqual('Pathogenic', calculate_interpretation('known'))
        self.assertEqual('Likely pathogenic', calculate_interpretation('likely'))
        self.assertEqual('Uncertain significance', calculate_interpretation('unknown'))
        self.assertEqual('other', calculate_interpretation('ambiguous'))
        self.assertEqual('', calculate_interpretation('fred'))

    def test_gather_attributes(self):
        copy_number = {
            '@equivocal': 'true',
            '@in-frame': 'unknown',
            '@other-gene': 'N/A',
            '@pos1': 'ch17:29557687-29887856',
            '@pos2': 'ch6:66426718-66427149',
            '@status': 'likely',
            '@supporting-read-pairs': 83,
            '@targeted-gene': 'NF1',
            '@type': 'truncation',
            'dna-evidence': {
                '@sample': 'SA-1612348'
            }
        }
        self.assertEqual(
            {'equivocal': 'true', 'supporting-read-pairs': 83},
            gather_attributes(copy_number))

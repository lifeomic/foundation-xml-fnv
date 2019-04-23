#!/usr/bin/env python
import argparse
import ast
import csv
import json
import logging

import xmltodict

logging.basicConfig(level=logging.INFO,
                    format='[%(asctime)s] %(levelname)s [%(name)s.%(funcName)s:%(lineno)d] %(message)s')
logger = logging.getLogger(__name__)


def read_xml(xml_file):
    with open(xml_file) as fd:
        return xmltodict.parse(fd.read())


def calculate_effect(fusion_type):
    if fusion_type == 'fusion':
        return 'fusion'
    if fusion_type == 'rearrangement':
        return 'rearrangement'
    if fusion_type == 'truncation':
        return 'truncation'
    if fusion_type == 'deletion':
        return 'deletion'
    if fusion_type == 'duplication':
        return 'duplication'
    if fusion_type == 'unknown':
        return 'unknown'

    logger.error('Failed to resolve fusion variant type: %s', fusion_type)
    return ''


def calculate_interpretation(status):
    if status == 'known':
        return 'Pathogenic'
    if status == 'likely':
        return 'Likely pathogenic'
    if status == 'unknown':
        return 'Uncertain significance'
    if status == 'ambiguous':
        return 'other'

    logger.error('Failed to resolve interpretation: %s', status)
    return ''


def gather_attributes(fusion):
    attributes = {}
    if '@type' in fusion.keys():
        attributes['type'] = fusion['@type']
    if '@equivocal' in fusion.keys():
        attributes['equivocal'] = fusion['@equivocal']
    if '@in-frame' in fusion.keys():
        attributes['in-frame'] = fusion['@in-frame']
    if '@supporting-read-pairs' in fusion.keys():
        attributes['supporting-read-pairs'] = fusion['@supporting-read-pairs']

    return attributes


def calculate_other_gene(other_gene):
    if other_gene == 'N/A':
        return ''
    return other_gene


def extract_fusion_variant(results_payload_dict):
    logger.info('Extracting fusion variants from xml')
    fusion_variant_list = {'FusionVariants': []}

    if 'rearrangements' in results_payload_dict['variant-report'].keys():
        if (results_payload_dict['variant-report']['rearrangements'] is not None and
                'rearrangement' in results_payload_dict['variant-report']['rearrangements'].keys()):
            variants_dict = results_payload_dict['variant-report']['rearrangements']['rearrangement']
            fusion_variants = variants_dict if isinstance(variants_dict, list) else [variants_dict]

            for fusion_variant in fusion_variants:
                fusion_variant_value = {'sample_id': fusion_variant['dna-evidence']['@sample'],
                                        'gene1': fusion_variant['@target-gene'],
                                        'gene2': calculate_other_gene(fusion_variant['@other-gene']),
                                        'effect': calculate_effect(fusion_variant['@type']),
                                        'chromosome1': fusion_variant['@pos1'].split(":")[0],
                                        'start_position1': int(fusion_variant['@pos1'].split(":")[1].split("-")[0]),
                                        'end_position1': int(fusion_variant['@pos1'].split(":")[1].split("-")[1]),
                                        'chromosome2': fusion_variant['@pos2'].split(":")[0],
                                        'start_position2': int(fusion_variant['@pos2'].split(":")[1].split("-")[0]),
                                        'end_position2': int(fusion_variant['@pos2'].split(":")[1].split("-")[1]),
                                        'interpretation': calculate_interpretation(fusion_variant['@status']),
                                        'sequence_type': 'somatic',
                                        'attributes': gather_attributes(fusion_variant)}
                fusion_variant_list['FusionVariants'].append(ast.literal_eval(json.dumps(fusion_variant_value)))

    return fusion_variant_list


def write_fusions_to_fnv(fnv_dict, args):
    logger.info('Saving fusion variants to fnv file')

    with open(args.out_file, 'w') as csvfile:
        csv_writer = csv.writer(csvfile, delimiter=',')
        csv_writer.writerow(['sample_id', 'gene1', 'gene2', 'effect', 'chromosome1', 'start_position1', 'end_position1',
                             'chromosome2', 'start_position2', 'end_position2', 'interpretation',
                             'sequence_type', 'attributes'])
        for fnv in fnv_dict['FusionVariants']:
            csv_writer.writerow([fnv['sample_id'], fnv['gene1'], fnv['gene2'],
                                 fnv['effect'], fnv['chromosome1'], fnv['start_position1'],
                                 fnv['end_position1'], fnv['chromosome2'], fnv['start_position2'],
                                 fnv['end_position2'], fnv['interpretation'],
                                 fnv['sequence_type'], fnv['attributes']])


def main():
    parser = argparse.ArgumentParser(
        prog='foundation-xml-fnv',
        description='Extracts fusion variants information from FoundationOne XML reports into CSV resources.')
    parser.add_argument('-x, --xml', dest='xml_file',
                        required=True, help='Path to the XML file')
    parser.add_argument('-o, --output', dest='out_file',
                        required=True, help='Path to write the FNV file')
    args = parser.parse_args()

    logger.info('Extracting fusion variants from XML and into FNV with args: %s', json.dumps(args.__dict__))
    xml_dict = read_xml(args.xml_file)

    fnv_dict = extract_fusion_variant(
        xml_dict['rr:ResultsReport']['rr:ResultsPayload'])

    write_fusions_to_fnv(fnv_dict, args)


if __name__ == '__main__':
    main()

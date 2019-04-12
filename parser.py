import logging
import os
import datetime
import time

FILE_NOT_FOUND_ERROR = 'Cannot find input file: {}'  # error message constant

# configure logger
logging.basicConfig(format='%(asctime)s %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s', level=logging.INFO)
_logger = logging.getLogger('3DTS-parser')

# change following parameters accordingly
SOURCE_NAME = '3dts'  # source name that appears in the api response
FILENAME = '3dts_loci_scores.sorted.tsv'  # name of the file to read
DELIMITER = '\t'  # the delimiter that separates each field


def version(self):
    return 'v5'


def _inspect_file(filename: str) -> int:
    _logger.info('start counting file lines')
    i = 0
    with open(filename) as f:
        for i, l in enumerate(f):
            pass
    _logger.info(f'total lines: {i+1}')
    return i + 1


def load_data(data_folder: str):
    """
    Load data from a specified file path. Parse each line into a dictionary according to the schema.
    Then process each dict by normalizing data format, remove null fields (optional).
    Append each dict into final result using its id.

    :param data_folder: the path(folder) where the data file is stored
    :return: a generator that yields data.
    """
    input_file = os.path.join(data_folder, FILENAME)
    # raise an error if file not found
    if not os.path.exists(input_file):
        _logger.error(FILE_NOT_FOUND_ERROR.format(input_file))
        raise FileExistsError(FILE_NOT_FOUND_ERROR.format(input_file))

    file_lines = _inspect_file(input_file)  # get total lines so that we can indicate progress in next step

    with open(input_file, 'r') as file:
        _logger.info(f'start reading file: {FILENAME}')
        count = 0
        skipped = []
        start_time = time.time()
        for line in file:
            count += 1
            ratio = count / file_lines
            time_left = datetime.timedelta(seconds=(time.time() - start_time) * (1 - ratio) / ratio)
            # format to use 2 decimals for progress
            if count % 10000 == 0:   # show progress every 500k records
                _logger.info(f'reading line {count} ({(ratio * 100):.2f}%), #skipped {len(skipped)}, estimated time left: {time_left}')

            if line.startswith('#') or line.strip() == '':
                skipped.append(line)
                continue  # skip commented/empty lines

            try:
                (chrom, start, end, score, pdb_id, pdb_chain, uniprot_feature_name, pdb_residue_min,
                 pdb_residue_max) = line.strip().split(DELIMITER)  # unpack according to schema
            except ValueError:
                _logger.error(f'failed to unpack line {count}: {line}')
                _logger.error(f'got: {line.strip().split(DELIMITER)}')
                skipped.append(line)
                continue  # skip error line

            try:  # parse each field if necessary (format, enforce datatype etc.)
                chrom = chrom.replace('chr', '')
                start = int(start)
                end = int(end)
                score = float(score)
            except ValueError as e:
                _logger.error(f'failed to cast type for line {count}: {e}')
                skipped.append(line)
                continue  # skip error line

            _id = f'chr{chrom}:g.{start}_{end}'  # define id

            variant = {
                'chrom': chrom,
                'start': start,
                'end': end,
                'scores': [
                    {
                        'score': score,
                        'pdb_id': pdb_id,
                        'pdb_chain': pdb_chain,
                        'uniprot_feature_name': uniprot_feature_name,
                        'pdb_residue_min': pdb_residue_min,
                        'pdb_residue_max': pdb_residue_max,
                    },
                ]
            }

            yield {  # commit an entry by yielding
                "_id": _id,
                SOURCE_NAME: variant
            }
        _logger.info(f'parse completed, {len(skipped)}/{file_lines} lines skipped.')
        for x in skipped:
            _logger.info(f'skipped line: {x.strip()}')

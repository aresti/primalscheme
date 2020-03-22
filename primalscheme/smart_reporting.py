import logging
import os
import pickle
from Bio import SeqIO, SeqUtils
from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import FeatureLocation, SeqFeature
from reportlab.lib import colors
from .smart import SMARTplexScheme

logger = logging.getLogger('Primal Log')

class SMARTplexReporter(SMARTplexScheme):
    """Reporting methods to extend MultiplexScheme"""

    def write_bed(self, path='./'):
        logger.info('Writing BED')
        filepath = os.path.join(path, '{}.bed'.format(self.prefix))
        with open(filepath, 'w') as bedhandle:
            for r in self.regions:
                pass

    def write_tsv(self, path='./'):
        logger.info('Writing TSV')
        filepath = os.path.join(path, '{}.tsv'.format(self.prefix))
        with open(filepath, 'w') as tsvhandle:
            print(*['name', 'pool', 'seq', 'length', '%gc', 'tm (use 65)'], sep='\t', file=tsvhandle)
            for r in self.regions:
                pass

    def write_pickle(self, path='./'):
        logger.info('Writing pickles')
        filepath = os.path.join(path, '{}.pickle'.format(self.prefix))
        with open(filepath, 'wb') as pickleobj:
            pickle.dump(self.regions, pickleobj)

    def write_refs(self, path='./'):
        logger.info('Writing references')
        filepath = os.path.join(path, '{}.fasta'.format(self.prefix))
        with open(filepath, 'w') as refhandle:
            SeqIO.write(self.references, filepath, 'fasta')

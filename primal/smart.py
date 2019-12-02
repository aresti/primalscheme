import logging
import primer3
from primal import settings
from .exceptions import MaxGapReached, NoSuitablePrimers
from .models import CandidatePrimer, CandidatePrimerPair, Region

logger = logging.getLogger('Primal Log')

class SMARTplexScheme(object):
    """A complete SMART-plex primer scheme."""

    def __init__(self, references, amplicon_length, max_candidates, prefix='PRIMAL_SCHEME'):
        self.references = references
        self.amplicon_length = amplicon_length
        self.max_candidates = max_candidates
        self.prefix = prefix
        self.regions = []

        self.run()

    @property
    def primary_reference(self):
        return self.references[0]

    def run(self):
        regions = []
        region_num = 0

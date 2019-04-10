class NoSuitablePrimers(Exception):
    """No suitable primer found."""

    def __init__(self, message):
        self.message = message


class MaxGapReached(Exception):
    """Maximum gap exceeded, increase amplicon length."""

    def __init__(self, message):
        self.message = message

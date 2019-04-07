class Error(Exception):
    """Base exception class"""

    pass


class NoSuitableError(Error):
    """No suitable primer found."""

    def __init__(self, message):
        self.message = message


class MaxGapError(Exception):
    """Maximum gap exceeded, increase amplicon length."""

    def __init__(self, message):
        self.message = message

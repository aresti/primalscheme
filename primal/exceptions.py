class NoSuitableException(Exception):
    """No suitable primer found."""

    pass


class MaxGapException(Exception):
    """Maximum gap exceeded, increase amplicon length."""

    pass

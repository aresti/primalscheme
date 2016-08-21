from __future__ import absolute_import
# Copyright (c) 2010-2016 openpyxl

import sys

VER = sys.version_info

from math import isnan
from .numbers import NUMERIC_TYPES

if VER[0] == 3:
    basestring = str
    unicode = str
    from io import BufferedReader
    file = BufferedReader
    from io import BufferedRandom
    tempfile = BufferedRandom
    bytes = bytes
else:
    basestring = basestring
    unicode = unicode
    file = file
    tempfile = file
    bytes = str


def safe_string(value):
    """Safely and consistently format numeric values"""
    if isinstance(value, NUMERIC_TYPES):
        if isnan(value):
            value = ""
        else:
            value = "%.16g" % value
    elif value is None:
        value = "none"
    elif not isinstance(value, basestring):
        value = str(value)
    return value

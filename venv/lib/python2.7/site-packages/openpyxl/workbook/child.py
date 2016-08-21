from __future__ import absolute_import
# Copyright (c) 2010-2016 openpyxl

import re
from openpyxl.compat import unicode

"""
Base class for worksheets, chartsheets, etc. that can be added to workbooks
"""

INVALID_TITLE_REGEX = re.compile(r'[\\*?:/\[\]]')


def avoid_duplicate_name(names, value):
    """
    Naive check to see whether name already exists.
    If name does exist suggest a name using an incrementer
    """
    if value in names:
        names = ",".join(names)
        sheet_title_regex = re.compile("(?P<title>%s)(?P<count>\d*),?" % re.escape(value))
        matches = sheet_title_regex.findall(names)
        if matches:
            # use name, but append with the next highest integer
            counts = [int(idx) for (t, idx) in matches if idx.isdigit()]
            highest = 0
            if counts:
                highest = max(counts)
            value = "{0}{1}".format(value, highest + 1)
    return value


class _WorkbookChild(object):

    __title = ""
    __parent = None
    _default_title = "Sheet"

    def __init__(self, parent=None, title=None):
        self.__parent = parent
        self.title = title or self._default_title


    def __repr__(self):
        return u'<{0} "{1}">'.format(self.__class__.__name__, self.title)


    @property
    def parent(self):
        return self.__parent


    @property
    def encoding(self):
        return self.__parent.encoding


    @property
    def title(self):
        return self.__title


    @title.setter
    def title(self, value):
        """
        Set a sheet title, ensuring it is valid.
        Limited to 31 characters, no special characters.
        Duplicate titles will be incremented numerically
        """
        if not value:
            raise ValueError("Title must have at least one character")

        if hasattr(value, "decode"):
            if not isinstance(value, unicode):
                try:
                    value = value.decode("ascii")
                except UnicodeDecodeError:
                    raise ValueError("Worksheet titles must be unicode")

        m = INVALID_TITLE_REGEX.search(value)
        if m:
            msg = "Invalid character {0} found in sheet title".format(m.group(0))
            raise ValueError(msg)

        if self.title is not None and self.title != value:
            value = avoid_duplicate_name(self.parent.sheetnames, value)

        if len(value) > 31:
            raise ValueError('Maximum 31 characters allowed in sheet title')

        self.__title = value

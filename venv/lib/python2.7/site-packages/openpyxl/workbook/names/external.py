from __future__ import absolute_import
# Copyright (c) 2010-2016 openpyxl

import os

from openpyxl.descriptors import String, Strict
from openpyxl.packaging.relationship import Relationship
from openpyxl.xml.constants import (
    SHEET_MAIN_NS,
    REL_NS,
    PKG_REL_NS,
    EXTERNAL_LINK_NS,
    ARC_WORKBOOK,
)
from openpyxl.xml.functions import (
    fromstring,
    safe_iterator,
    Element,
    SubElement,
)


"""Manage links to external Workbooks"""


class ExternalBook(Strict):

    """
    Map the relationship of one workbook to another
    """

    Id = String()
    Type = "http://schemas.openxmlformats.org/officeDocument/2006/relationships/externalLinkPath"
    TargetMode = "External"
    Target = String()

    def __init__(self, Id, Target, TargetMode=None, Type=None):
        self.Id = Id
        self.Target = Target
        links = []

    def __iter__(self):
        for attr in ('Id', 'Type', 'TargetMode', 'Target'):
            value = getattr(self, attr)
            yield attr, value


class ExternalRange(Strict):

    """
    Map external named ranges
    NB. the specification for these is different to named ranges within a workbook
    See 18.14.5
    """

    name = String()
    refersTo = String(allow_none=True)
    sheetId = String(allow_none=True)

    def __init__(self, name, refersTo=None, sheetId=None):
        self.name = name
        self.refersTo = refersTo
        self.sheetId = sheetId


    def __iter__(self):
        for attr in ('name', 'refersTo', 'sheetId'):
            value = getattr(self, attr, None)
            if value is not None:
                yield attr, value


def parse_books(xml):
    tree = fromstring(xml)
    rels = tree.findall('{%s}Relationship' % PKG_REL_NS)
    for r in rels:
        return ExternalBook(**r.attrib)


def parse_ranges(xml):
    tree = fromstring(xml)
    book = tree.find('{%s}externalBook' % SHEET_MAIN_NS)
    if book is None:
        return
    names = book.find('{%s}definedNames' % SHEET_MAIN_NS)
    for n in safe_iterator(names, '{%s}definedName' % SHEET_MAIN_NS):
        yield ExternalRange(**n.attrib)


def detect_external_links(rels, archive):
    from openpyxl.reader.workbook import find_external_refs
    rels = dict(rels)

    for rId in find_external_refs(archive):
        rel = rels[rId]
        pth = os.path.split(rel['path'])
        f_name = pth[-1]
        dir_name = "/".join(pth[:-1])
        book_path = "{0}/_rels/{1}.rels".format (dir_name, f_name)
        book_xml = archive.read(book_path)
        Book = parse_books(book_xml)

        range_xml = archive.read(rel['path'])
        Book.links = list(parse_ranges(range_xml))
        yield Book


def write_external_link(links):
    """Serialise links to ranges in a single external worbook"""
    root = Element("{%s}externalLink" % SHEET_MAIN_NS)
    book =  SubElement(root, "{%s}externalBook" % SHEET_MAIN_NS, {'{%s}id' % REL_NS:'rId1'})
    external_ranges = SubElement(book, "{%s}definedNames" % SHEET_MAIN_NS)
    for l in links:
        external_ranges.append(Element("{%s}definedName" % SHEET_MAIN_NS, dict(l)))
    return root


def write_external_book_rel(book):
    """Serialise link to external file"""
    root = Element("Relationships", xmlns=PKG_REL_NS)
    rel = Relationship("", target=book.Target, targetMode=book.TargetMode, id="rId1")
    rel.type = book.Type
    root.append(rel.to_tree())
    return root

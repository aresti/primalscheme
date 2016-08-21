from __future__ import absolute_import
# Copyright (c) 2010-2016 openpyxl


import os.path

from openpyxl.comments import Comment
from openpyxl.xml.constants import (
    PACKAGE_WORKSHEET_RELS,
    COMMENTS_NS,
    PACKAGE_XL,
    )
from openpyxl.xml.functions import fromstring

from .properties import CommentSheet


def read_comments(ws, xml_source):
    """Given a worksheet and the XML of its comments file, assigns comments to cells"""
    root = fromstring(xml_source)
    comments = CommentSheet.from_tree(root)
    authors = comments.authors.author

    for comment in comments.commentList:
        author = authors[comment.authorId]
        ref = comment.ref
        comment = Comment(comment.content, author)

        ws.cell(coordinate=ref).comment = comment


def get_comments_file(worksheet_path, archive, valid_files):
    """Returns the XML filename in the archive which contains the comments for
    the spreadsheet with codename sheet_codename. Returns None if there is no
    such file"""
    sheet_codename = os.path.split(worksheet_path)[-1]
    rels_file = PACKAGE_WORKSHEET_RELS + '/' + sheet_codename + '.rels'
    if rels_file not in valid_files:
        return None
    rels_source = archive.read(rels_file)
    root = fromstring(rels_source)
    for i in root:
        if i.attrib['Type'] == COMMENTS_NS:
            comments_file = os.path.split(i.attrib['Target'])[-1]
            comments_file = PACKAGE_XL + '/' + comments_file
            if comments_file in valid_files:
                return comments_file
    return None

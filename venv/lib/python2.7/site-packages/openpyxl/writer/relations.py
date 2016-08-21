from __future__ import absolute_import
# Copyright (c) 2010-2016 openpyxl

from openpyxl.packaging.relationship import Relationship, RelationshipList


def write_rels(worksheet, comments_id=None):
    """Write relationships for the worksheet to xml."""

    rels = RelationshipList(worksheet._rels)

    # If there is an existing vml file that is preserved or extended then
    # create its relation.
    if worksheet.legacy_drawing is not None:
        rel = Relationship("vmlDrawing", id="anysvml", target='/' + worksheet.legacy_drawing)
        rels.append(rel)

    # Comments
    if worksheet._comment_count > 0:
        rel = Relationship(type="comments", id="comments",
                           target='/xl/comments%s.xml' % comments_id)
        rels.append(rel)

        if worksheet.legacy_drawing is None:
            rel = Relationship(type="vmlDrawing", id="anysvml",
                           target='/xl/drawings/commentsDrawing%s.vml' % comments_id)
            rels.append(rel)

    return rels.to_tree()

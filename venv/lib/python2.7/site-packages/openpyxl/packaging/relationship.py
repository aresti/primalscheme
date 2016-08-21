from __future__ import absolute_import
# Copyright (c) 2010-2016 openpyxl

import posixpath

from openpyxl.descriptors import (
    String,
    Set,
    NoneSet,
    Alias,
    Sequence,
)
from openpyxl.descriptors.serialisable import Serialisable

from openpyxl.xml.constants import REL_NS, PKG_REL_NS
from openpyxl.xml.functions import (
    Element,
    fromstring,
    tostring
)


class Relationship(Serialisable):
    """Represents many kinds of relationships."""
    # TODO: Use this object for workbook relationships as well as
    # worksheet relationships

    tagname = "Relationship"

    Type = String()
    type = Alias('Type')
    Target = String()
    target = Alias('Target')
    TargetMode = String(allow_none=True)
    targetMode = Alias('TargetMode')
    Id = String(allow_none=True)
    id = Alias('Id')


    def __init__(self,
                 type=None,
                 target=None,
                 targetMode=None,
                 id=None,
                 Id=None,
                 Type=None,
                 Target=None,
                 ):
        if type is not None:
            Type = "%s/%s" % (REL_NS, type)
        self.Type = Type
        if target is not None:
            Target = target
        self.Target = Target
        self.targetMode = targetMode
        if id is not None:
            Id = id
        self.Id = Id


class RelationshipList(Serialisable):

    tagname = "Relationships"

    Relationship = Sequence(expected_type=Relationship)


    def __init__(self, Relationship=()):
        self.Relationship = Relationship


    def append(self, value):
        values = self.Relationship[:]
        values.append(value)
        self.Relationship = values


    def __len__(self):
        return len(self.Relationship)


    def __bool__(self):
        return bool(self.Relationship)


    def __getitem__(self, key):
        for r in self.Relationship:
            if r.id == key:
                return r
        raise KeyError("Unknown relationship: {0}".format(key))


    def to_tree(self):
        tree = Element("Relationships", xmlns=PKG_REL_NS)
        for idx, rel in enumerate(self.Relationship, 1):
            if not rel.id:
                rel.id = "rId{0}".format(idx)
            tree.append(rel.to_tree())

        return tree


def get_dependents(archive, filename):
    """
    Normalise dependency file paths to absolute ones

    Relative paths are relative to parent object
    """
    src = archive.read(filename)
    node = fromstring(src)
    rels = RelationshipList.from_tree(node)
    folder = posixpath.dirname(filename)
    parent = posixpath.split(folder)[0]
    for r in rels.Relationship:
        if r.target.startswith("/"):
            r.target = r.target[1:]
            continue
        pth = posixpath.join(parent, r.target)
        r.target = posixpath.normpath(pth)
    return rels

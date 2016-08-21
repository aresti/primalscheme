from __future__ import absolute_import
# Copyright (c) 2010-2016 openpyxl

"""
File manifest
"""
import mimetypes
import os.path

from openpyxl.descriptors.serialisable import Serialisable
from openpyxl.descriptors import String, Sequence
from openpyxl.xml.functions import fromstring
from openpyxl.xml.constants import (
    ARC_CORE,
    ARC_CONTENT_TYPES,
    ARC_WORKBOOK,
    ARC_APP,
    ARC_THEME,
    ARC_STYLE,
    ARC_SHARED_STRINGS,
    EXTERNAL_LINK,
    THEME_TYPE,
    STYLES_TYPE,
    XLSX,
    XLSM,
    XLTM,
    XLTX,
    WORKSHEET_TYPE,
    COMMENTS_TYPE,
    SHARED_STRINGS,
    DRAWING_TYPE,
    CHART_TYPE,
    CHARTSHAPE_TYPE,
    CHARTSHEET_TYPE,
    CONTYPES_NS
)

# initialise mime-types
mimetypes.init()
mimetypes.add_type('application/xml', ".xml")
mimetypes.add_type('application/vnd.openxmlformats-package.relationships+xml', ".rels")
mimetypes.add_type("application/vnd.ms-office.activeX", ".bin")
mimetypes.add_type("application/vnd.openxmlformats-officedocument.vmlDrawing", ".vml")
mimetypes.add_type("image/x-emf", ".emf")


class FileExtension(Serialisable):

    tagname = "Default"

    Extension = String()
    ContentType = String()

    def __init__(self, Extension, ContentType):
        self.Extension = Extension
        self.ContentType = ContentType


    def __hash__(self):
        return hash((self.Extension, self.ContentType))


class Override(Serialisable):

    tagname = "Override"

    PartName = String()
    ContentType = String()

    def __init__(self, PartName, ContentType):
        self.PartName = PartName
        self.ContentType = ContentType


    def __hash__(self):
        return hash((self.PartName, self.ContentType))


DEFAULT_TYPES = [
    FileExtension("rels", "application/vnd.openxmlformats-package.relationships+xml"),
    FileExtension("xml", "application/xml"),
]

DEFAULT_OVERRIDE = [
    Override("/" + ARC_WORKBOOK, XLSX), # Workbook
    Override("/" + ARC_SHARED_STRINGS, SHARED_STRINGS), # Shared strings
    Override("/" + ARC_STYLE, STYLES_TYPE), # Styles
    Override("/" + ARC_THEME, THEME_TYPE), # Theme
    Override("/docProps/core.xml", "application/vnd.openxmlformats-package.core-properties+xml"),
    Override("/docProps/app.xml", "application/vnd.openxmlformats-officedocument.extended-properties+xml")
]


class Manifest(Serialisable):

    tagname = "Types"

    Default = Sequence(expected_type=FileExtension, unique=True)
    Override = Sequence(expected_type=Override, unique=True)

    __elements__ = ("Default", "Override")

    def __init__(self,
                 Default=(),
                 Override=(),
                 ):
        if not Default:
            Default = DEFAULT_TYPES
        self.Default = Default
        if not Override:
            Override = DEFAULT_OVERRIDE
        self.Override = Override


    @property
    def filenames(self):
        return [part.PartName for part in self.Override]


    @property
    def extensions(self):
        exts = set([os.path.splitext(part.PartName)[-1] for part in self.Override])
        return [(ext[1:], mimetypes.types_map[ext]) for ext in sorted(exts)]


    def to_tree(self):
        """
        Custom serialisation method to allow setting a default namespace
        """
        defaults = [t.Extension for t in self.Default]
        for ext, mime in self.extensions:
            if ext not in defaults:
                mime = FileExtension(ext, mime)
                self.Default.append(mime)
        tree = super(Manifest, self).to_tree()
        tree.set("xmlns", CONTYPES_NS)
        return tree


def write_content_types(workbook, as_template=False, exts=None):

    manifest = Manifest()

    if exts is not None:
        for ext in exts:
            ext = os.path.splitext(ext)[-1]
            mime = mimetypes.types_map[ext]
            fe = FileExtension(ext[1:], mime)
            if fe not in manifest.Default:
                manifest.Default.append(fe)

    if workbook.vba_archive:
        node = fromstring(workbook.vba_archive.read(ARC_CONTENT_TYPES))
        manifest = Manifest.from_tree(node)
        del node
        partnames = [t.PartName for t in manifest.Override]
        for override in DEFAULT_OVERRIDE:
            if override.PartName not in partnames:
                manifest.Override.append(override)

    # templates
    for part in manifest.Override:
        if part.PartName == "/" + ARC_WORKBOOK:
            ct = as_template and XLTX or XLSX
            if workbook.vba_archive:
                ct = as_template and XLTM or XLSM
            part.ContentType = ct


    drawing_id = 0
    chart_id = 0
    comments_id = 0

    # ugh! can't we get this from the zip archive?
    # worksheets
    for sheet_id, sheet in enumerate(workbook.worksheets):
        name = '/xl/worksheets/sheet%d.xml' % (sheet_id + 1)
        manifest.Override.append(Override(name, WORKSHEET_TYPE))

        if sheet._charts or sheet._images:
            drawing_id += 1
            name = '/xl/drawings/drawing%d.xml' % drawing_id
            manifest.Override.append(Override(name, DRAWING_TYPE))


            for chart in sheet._charts:
                chart_id += 1
                name = '/xl/charts/chart%d.xml' % chart_id
                manifest.Override.append(Override(name, CHART_TYPE))

        if sheet._comment_count > 0:
            comments_id += 1
            vml = FileExtension("vml", mimetypes.types_map[".vml"])
            if vml not in manifest.Default:
                manifest.Default.append(vml)
            name = '/xl/comments%d.xml' % comments_id
            manifest.Override.append(Override(name, COMMENTS_TYPE))


    # chartsheets
    for sheet_id, sheet in enumerate(workbook.chartsheets, sheet_id+1):
        name = '/xl/chartsheets/sheet%d.xml' % (sheet_id)
        manifest.Override.append(Override(name, CHARTSHEET_TYPE))

        if sheet._charts:
            drawing_id += 1
            name = '/xl/drawings/drawing%d.xml' % drawing_id
            manifest.Override.append(Override(name, DRAWING_TYPE))

            for chart in sheet._charts:
                chart_id += 1
                name = '/xl/charts/chart%d.xml' % chart_id
                manifest.Override.append(Override(name, CHART_TYPE))

    #external links
    for idx, _ in enumerate(workbook._external_links, 1):
        name = '/xl/externalLinks/externalLink{0}.xml'.format(idx)
        manifest.Override.append(Override(name, EXTERNAL_LINK))

    return manifest

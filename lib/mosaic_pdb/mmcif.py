# -*- coding: utf-8 -*-

"""mmCIF parser

The functions in this module parse mmCIF files and return the contained
data items. There are four different routines depending on the level of
access that the application requires.

.. moduleauthor:: Konrad Hinsen <konrad.hinsen@cnrs-orleans.fr>

"""

#-----------------------------------------------------------------------------
#       Copyright (C) 2013 The Mosaic Development Team
#
#  Distributed under the terms of the BSD License.  The full license is in
#  the file LICENSE.txt, distributed as part of this software.
#-----------------------------------------------------------------------------

# The first stage of this parser (tokenizer) was inspired by the mmCIF parser
# in PyMMLIB (http://pymmlib.sourceforge.net/). The regular expression for
# identifying tokens was taken directly from there. The higher-level part of
# the parser is completely different, however. The PyMMLIB parser loads
# all the file contents into memory. This turned out to consume too much
# memory for large structure factor file. The parser in this module stores
# nothing at all. It is written as an iterator yielding mmCIF data items
# to the calling routine. The highest-level interface provides a filter
# for selecting data items and feeds tables to the calling program line
# by line through a callback object. This ensures that only the required
# data is stored, and directly in the form that the application will use.

import re
import gzip
import os
import sys

import numpy as np

# Python 2/3 compatibility issues
if sys.version_info[0] == 2:
    from urllib import urlretrieve
    bytes2text = lambda x: x
else:
    from urllib.request import urlretrieve
    from io import TextIOWrapper
    def bytes2text(stream):
        return TextIOWrapper(stream, encoding="utf8")

# Identifiers used in the return values of the parsers.
TOKEN = 'token'
DATA_LABEL = 'data_label'
DATA_VALUE = 'data_value'
KEYWORD = 'keyword'

LABEL_OR_KEYWORD = 'label_or_keyword'
VALUE = 'value'
LOOP_LABELS = 'loop_labels'
LOOP_VALUES = 'loop_values'

DATA = 'data'
TABLE_HEADER = 'table_header'
TABLE_DATA = 'table_data'

# A regular expression that encodes mmCIF syntax.
_token_regexp = re.compile(
    r"(?:"
     "(?:_(.+?)[.](\S+))"               "|"  # _section.subsection
     "(?:['\"](.*?)(?:['\"]\s|['\"]$))" "|"  # quoted strings
     "(?:\s*#.*$)"                      "|"  # comments
     "(\S+)"                                 # unquoted tokens
     ")")


# mmCIF-specific Exception classes

class MMCIFError(Exception):
    pass


class MMCIFSyntaxError(MMCIFError):

    def __init__(self, text, line_number=None):
        if line_number is None:
            # This is a workaround for the __setstate__ method in Exception
            # that effectively requires exceptions to accept a single
            # argument, otherwise a bug occurs during unpickling.
            MMCIFError.__init__(self, text)
            return
        self.line_number = line_number
        MMCIFError.__init__(self, "Line %d: %s" % (line_number, text))


#
# The parser object. It provides interfaces at different levels.
#
class MMCIFParser(object):

    """
    Parser for mmCIF files

    The file to be parsed is specified when the parser object is created.
    One of the parser functions/generators may then be called to access
    the file contents.
    """

    def __init__(self, file_name=None, file_object=None, pdb_code=None):
        """
        Specify the mmCIF file to be loaded. Only one of the four
        keyword parameters may be given a value.

        :param file_name: the name of a file
        :type file_name: str
        :param file_object: a file object
        :type file_object: file
        :param pdb_code: the PDB code for a structure file, which is
            taken from a public or local PDB repository
        :type pdb_code: str
        """
        self.line_number = 0
        if file_name is not None:
            assert file_object is None
            assert pdb_code is None
            self.file_object = open(file_name)
        elif file_object is not None:
            assert pdb_code is None
            self.file_object = file_object
        elif pdb_code is not None:
            self.file_object = mmcif_files.getFile(pdb_code)
        else:
            raise ValueError("No input file given")

    def parseLowLevel(self):
        """
        An iterator that yields the contents of the mmCIF file in the
        form of (type, data) pairs. The type can be KEYWORD, DATA_LABEL,
        or DATA_VALUE.
        """
        file_iter = iter(self.file_object)
        while True:
            line = next(file_iter)
            self.line_number += 1

            ## skip comments
            if line.startswith("#"):
                continue

            ## semi-colon multi-line strings
            if line.startswith(";"):
                lmerge = [line[1:]]
                while True:
                    line = next(file_iter)
                    self.line_number += 1
                    if line.startswith(";"):
                        break
                    lmerge.append(line)

                lmerge[-1] = lmerge[-1].rstrip()
                yield (DATA_VALUE, "".join(lmerge))
                continue

            ## split line into tokens
            for match in _token_regexp.finditer(line):
                label1, label2, string, token = match.groups()
                if label1 is not None and label2 is not None:
                    yield DATA_LABEL, (label1, label2)
                elif string is not None:
                    yield DATA_VALUE, string
                elif token is not None:
                    token_parts = token.split('_')
                    if len(token_parts) == 1 \
                       or token_parts[0].lower() not in ("data", "loop",
                                                         "global", "save",
                                                         "stop"):
                        yield DATA_VALUE, token
                    else:
                        yield KEYWORD, (token_parts[0].lower(),
                                        '_'.join(token_parts[1:]))

    def parse(self):
        """
        An iterator that yields the contents of the mmCIF file in the
        form of (type, data) pairs. The type can be KEYWORD (data
        is an optional label), DATA (data is a triple of (category label,
        item label, value)), TABLE_HEADER (data is a list of the category and
        item labels in the table) or TABLE_DATA (data is a list containing the
        data items for one row of the table).
        """
        iterator = self.parseLowLevel()
        while True:
            item_type, item = next(iterator)
            if item_type is KEYWORD and item[0] == "data":
                yield item_type, item
                break

        state = LABEL_OR_KEYWORD
        while True:
            item_type, item = next(iterator)

            if state is LABEL_OR_KEYWORD:
                if item_type is DATA_LABEL:
                    label1, label2 = item
                    state = VALUE
                elif item_type is KEYWORD:
                    if item[0] == "loop":
                        loop_labels = []
                        state = LOOP_LABELS
                    else:
                        yield item_type, item
                else:
                    raise MMCIFSyntaxError("Expected data label or keyword",
                                           self.line_number)

            elif state is VALUE:
                if item_type is DATA_VALUE:
                    state = LABEL_OR_KEYWORD
                    yield DATA, (label1, label2, item)
                else:
                    raise MMCIFSyntaxError("Expected data value "
                                           "for label %s.%s"
                                           % (label1, label2),
                                           self.line_number)

            elif state is LOOP_LABELS:
                if item_type is DATA_LABEL:
                    if loop_labels and loop_labels[0][0] != item[0]:
                        # The label does not belong to the loop category.
                        # meaning that the loop is empty and terminated.
                        label1, label2 = item
                        state = VALUE
                    else:
                        loop_labels.append(item)
                elif item_type is DATA_VALUE:
                    loop_data = [item]
                    state = LOOP_VALUES
                    yield TABLE_HEADER, loop_labels
                    if len(loop_labels) == 1:
                        yield TABLE_DATA, loop_data
                        loop_data = []
                else:
                    raise MMCIFSyntaxError("Expected label or value in loop",
                                           self.line_number)

            elif state is LOOP_VALUES:
                if item_type is DATA_VALUE:
                    loop_data.append(item)
                    if len(loop_data) == len(loop_labels):
                        yield TABLE_DATA, loop_data
                        loop_data = []
                else:
                    if len(loop_data) > 0:
                        raise MMCIFSyntaxError("Extraneous data in loop:" +
                                               str(loop_data),
                                               self.line_number)
                    if item_type is DATA_LABEL:
                        label1, label2 = item
                        state = VALUE
                    elif item_type is KEYWORD:
                        if item[0] == "loop":
                            loop_labels = []
                            state = LOOP_LABELS
                        else:
                            yield item_type, item
                    else:
                        raise MMCIFSyntaxError("Expected data label or loop",
                                               self.line_number)

    def parseAndSelect(self, categories, data=0):
        """
        An iterator that yields a subset of the contents of the mmCIF file
        in the form of (type, data) pairs. The return values are the same
        as for parse. However, only items corresponding to the selection
        are yielded. The selection consists of a list of categories and of
        the specification of a data set by number (0 = first) or name.
        Note that most mmCIF files contain only a single data set.

        :param categories: a sequence of category names
        :type categories: sequence of str
        :param data: the selected data set, either by number (0 being the
            first data set) or by name
        :type data: int or str
        """
        dataset = -1
        dataset_name = None
        for item_type, item in self.parse():

            if item_type is KEYWORD:
                if item[0] == 'data':
                    dataset += 1
                    dataset_name = item[1]
                    return_data = data == dataset or data == dataset_name
                else:
                    raise MMCIFError("Keyword %s not yet implemented"
                                     % item[0])

            elif item_type is DATA:
                if item[0] in categories and return_data:
                    yield DATA, item

            elif item_type is TABLE_HEADER:
                keep_table = False
                for label1, label2 in item:
                    if label1 in categories:
                        keep_table = True
                if keep_table and return_data:
                    yield TABLE_HEADER, item

            elif item_type is TABLE_DATA:
                if keep_table and return_data:
                    yield TABLE_DATA, item

            else:
                raise MMCIFSyntaxError("Unexpected item type %s"
                                       % str(item_type),
                                       self.line_number)

    def parseToObjects(self, data=0, **categories):
        """
        Parse the file and store the selected data using the provided
        keyword parameters. Each keyword argument specifies one category
        of data items to be processed; categories for which there is no
        keyword parameter are ignored. The value of the keyword argument
        should be a dictionary for standard data items (which are stored
        in the dictionary) and a callable object for tables. The object
        is called once for each row in the table. The two arguments
        given are 1) a dictionary mapping item labels to indices into
        the data list and 2) the items in the row as a list.

        :param data: the selected data set, either by number (0 being the
            first data set) or by name
        :type data: int or str
        :param categories: the category-specific handlers
        :type categories: dict
        """
        accumulator = None
        for item_type, item in self.parseAndSelect(categories.keys(), data):

            if item_type is DATA:
                label1, label2, value = item
                handler = categories[label1]
                if isinstance(handler, dict):
                    if accumulator is not None:
                        accumulator.deliver()
                        accumulator = None
                    handler[label2] = value
                else:
                    if accumulator is not None \
                           and not accumulator.sameLabel(label1):
                        accumulator.deliver()
                        accumulator = None
                    if accumulator is None:
                        accumulator = DataAccumulator(label1, handler)
                    accumulator.add(label2, value)

            else:

                if accumulator is not None:
                    accumulator.deliver()
                    accumulator = None

                if item_type is TABLE_HEADER:
                    indices = {}
                    for i, (label1, label2) in enumerate(item):
                        indices[label2] = i
                    handler = categories[label1]

                elif item_type is TABLE_DATA:
                    handler(indices, item)

                else:
                    raise MMCIFSyntaxError("Unexpected item type %s"
                                           % str(item_type),
                                           self.line_number)


class DataAccumulator(object):

    def __init__(self, label, handler):
        self.label = label
        self.handler = handler
        self.fields = []
        self.values = []

    def add(self, field, value):
        self.fields.append(field)
        self.values.append(value)

    def deliver(self):
        self.handler(dict((f, i) for i, f in  enumerate(self.fields)),
                     self.values)

    def sameLabel(self, label):
        return self.label == label


#
# Repository of PDB files
#

class PDBFileCollection(object):

    """
    A PDBFileCollection describes a database directory in a local PDB
    repository, containing either PDB files, mmCIF files, or structure
    factor files.
    """

    def __init__(self, base_path, filename_pattern, url_pattern):
        """
        :param base_path: the path to the database directory, i.e. the
                          path containing the two-letter subdirectories
        :type base_path: str
        :param filename_pattern: the pattern of the file names in the
                                 database, containing %s at the place of
                                 the four-letter PDB code.
        :type filename_pattern: str
        :param url_pattern: the pattern of the URLs for download from the
                            PDB server, containing %s at the place of
                            the four-letter PDB code.
        :type url_pattern: str
        """
        self.base_path = base_path
        self.is_local = base_path is not None
        self.filename_pattern = filename_pattern
        self.url_pattern = url_pattern
        self.pdb_code_index = self.filename_pattern.find('%s')

    def getFilename(self, pdb_code):
        """
        :param pdb_code: the four-letter PDB code
        :type pdb_code: str
        :return: the corresponding file name
        :rtype: str
        """
        assert len(pdb_code) == 4, "Invalid PDB code " + repr(pdb_code)
        if self.base_path is None:
            raise IOError("Directory path undefined")
        pdb_code = pdb_code.upper()
        subdir = pdb_code[1:3]
        return os.path.join(self.base_path, subdir,
                            self.filename_pattern % pdb_code)

    def fileExists(self, pdb_code):
        """
        :param pdb_code: the four-letter PDB code
        :type pdb_code: str
        :return: True if there is a corresponding file
        :rtype: bool
        """
        return os.path.exists(self.getFilename(pdb_code))

    def getFile(self, pdb_code):
        """
        :param pdb_code: the four-letter PDB code
        :type pdb_code: str
        :return: the corresponding file
        :rtype: file
        """
        if not self.is_local:
            assert len(pdb_code) == 4, "Invalid PDB code " + repr(pdb_code)
            pdb_code = pdb_code.upper()
            if self.url_pattern is None:
                raise IOError("No URL pattern for PDB repository")
            url = self.url_pattern % pdb_code
            filename, headers = urlretrieve(url)
            if url.endswith('.gz'):
                return bytes2text(gzip.GzipFile(filename))
            else:
                return open(filename)
        filename = self.getFilename(pdb_code)
        if filename.endswith('.gz'):
            return bytes2text(gzip.GzipFile(filename))
        else:
            return open(filename)

    def __iter__(self):
        """
        :return: a generator yielding the PDB codes for all the files
                 in the collection
        :rtype: generator
        """
        for dirpath, dirname, filenames in os.walk(self.base_path):
            for filename in filenames:
                pdb_code = filename[self.pdb_code_index:
                                    self.pdb_code_index + 4]
                if self.filename_pattern % pdb_code == filename:
                    yield pdb_code


pdb_mmcif_path = os.environ.get('PDB_MMCIF_PATH', None)
mmcif_files = PDBFileCollection(pdb_mmcif_path, '%s.cif.gz',
                                'http://www.rcsb.org/pdb/files/%s.cif.gz')

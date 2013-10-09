# -*- coding: utf-8 -*-

#-----------------------------------------------------------------------------
#       Copyright (C) 2013 The Mosaic Development Team
#
#  Distributed under the terms of the BSD License.  The full license is in
#  the file LICENSE.txt, distributed as part of this software.
#-----------------------------------------------------------------------------

import string
import sys
import weakref

import numpy as N


# Python 2/3 compatibility issues
if sys.version_info[0] == 2:

    def ascii(string):
        return string

    def py_str(byte_string):
        if isinstance(byte_string, N.ndarray):
            return str(byte_string)
        else:
            assert isinstance(byte_string, str)
            return byte_string

    def isstring(s):
        return isinstance(s, basestring)

    xml_encoding = "utf-8"

else:

    def ascii(string):
        return bytes(string, 'ASCII')

    def py_str(byte_string):
        if isinstance(byte_string, N.ndarray):
            byte_string = bytes(byte_string)
        assert isinstance(byte_string, bytes)
        return byte_string.decode('ASCII')

    def isstring(s):
        return isinstance(s, str)

    xml_encoding = "unicode"


# Test for ASCII string without control characters.
# This doesn't look very efficient, but it works under
# both Python 2 and Python 3. Under Python 2, it refuses
# objects of type unicode.
ascii_chars = string.ascii_letters + string.digits + string.punctuation + ' '


def isascii(s):
    return isinstance(s, str) and all(c in ascii_chars for c in s)


def uint_for_max_value(m):
    if m < 1 << 8:
        return N.uint8
    elif m < 1 << 16:
        return N.uint16
    elif m < 1 << 32:
        return N.uint32
    else:
        return N.uint64


class SymbolDict(dict):

    def __getitem__(self, item):
        i = dict.get(self, item, None)
        if i is None:
            i = len(self)
            self[item] = i
        return i


class MethodRegister(object):

    def __init__(self):
        self.register_cls = {}
        self.register_str = {}

    def __call__(self, type_id):
        def register(method):
            if isinstance(type_id, str):
                self.register_str[type_id] = method
            else:
                self.register_cls[type_id] = method
            return method
        return register

    def __getitem__(self, type_id):
        if isinstance(type_id, str):
            return self.register_str[type_id]
        else:
            for t, m in self.register_cls.items():
                if issubclass(type_id, t):
                    return m
            return None


class AbstractFactory(object):

    def __init__(self):
        self.memo = weakref.WeakKeyDictionary()

    def __call__(self, obj):
        try:
            return self.memo[obj]
        except KeyError:
            h = self.handler[type(obj)]
            if h is None:
                raise TypeError("No factory routine for %s"
                                % str(type(obj)))
            new_obj = h(self, obj)
            self.memo[obj] = new_obj
            return new_obj

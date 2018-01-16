# -*- coding: utf-8 -*-

import pathlib


def make_relpath(args):
    path = pathlib.Path(__file__).parent
    for arg in args:
        path /= arg
    return path

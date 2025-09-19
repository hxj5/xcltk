# base.py - basic utils.

import os

def assert_e(path):
    assert path is not None and os.path.exists(path)

def assert_n(x):
    assert x is not None and len(x) > 0

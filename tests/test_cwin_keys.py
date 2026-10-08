"""
input keywords of seedname.cwin must be documented consistently.

CWin._default() defines the accepted keywords; they are also listed in
util/header.py (shown by `pw2cw -i`), the CWin.read docstring and docs/src/format/cwin.md.
"""

import os
import re

import pytest

from symclosestwannier.cw.cwin import CWin
from symclosestwannier.util.header import cwin_info

from conftest import TOPDIR

KEYS = set(CWin._default().keys())


# ==================================================
def listed_keys(text):
    return set(re.findall(r"^\s*- (\w+)\s*:", text, re.M))


# ==================================================
@pytest.mark.parametrize(
    "name, keys",
    [
        ("util/header.py (pw2cw -i)", set(cwin_info.keys())),
        ("CWin.read docstring", listed_keys(CWin.read.__doc__)),
        ("docs/src/format/cwin.md", listed_keys(open(os.path.join(TOPDIR, "docs", "src", "format", "cwin.md")).read())),
    ],
)
def test_keys_are_documented(name, keys):
    assert sorted(KEYS - keys) == [], f"keywords missing in {name}"
    assert sorted(keys - KEYS) == [], f"unknown keywords in {name}"


# ==================================================
def test_header_has_no_duplicated_keys():
    with open(os.path.join(TOPDIR, "symclosestwannier", "util", "header.py")) as fp:
        text = fp.read().split("\nwin_info = {")[0]

    keys = re.findall(r'^\s*"(\w+)":', text, re.M)
    duplicated = sorted({k for k in keys if keys.count(k) > 1})

    assert duplicated == []

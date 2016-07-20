"""
fast5.__init__.py
(c) 2016: Matei David, Ontario Institute for Cancer Research
MIT License
"""

from .version import __version__
from fast5 import *

__version_info__ = tuple([int(num) for num in __version__.split('.')]) 

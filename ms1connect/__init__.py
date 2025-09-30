"""
MS1Connect: A tool that scores the similarity between a pair of mass spectrometry runs.

MS1Connect solves the challenging problem of comparing mass spectrometry data acquired 
under different experimental protocols by framing it as a maximum bipartite matching 
problem and using only data from intact peptide (MS1) scans.
"""

__version__ = "0.1.0"
__author__ = "Andy"
__email__ = "22140743+bmx8177@users.noreply.github.com"

# Import main functionality
from .core import *
from .ms1_feature_detection import *
from .edge_matrix import *
from .pairwise_edge_matrix import *
from .plots import *

__all__ = [
    '__version__',
    '__author__', 
    '__email__'
]

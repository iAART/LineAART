import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from aart_func import *
from params import * 

obsint.Redshifts(spin_case,i_case,sub_kep,betar,betaphi)
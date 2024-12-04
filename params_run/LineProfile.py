import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from aart_func import *
from params import * 

lp.LineProfile(spin_case,i_case,sub_kep,betar,betaphi,sigma,r_0,E_s,r_min,r_max,E_binwidth)

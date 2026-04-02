# Reactor 2D 
from fipy import *
import numpy as np
import matplotlib.pyplot as plt
from fipy.tools import numerix
import pandas as pd
from fipy.terms import ImplicitSourceTerm
import hvplot.pandas
import panel as pn
import pickle


#Random equations for non-NaN results
def safe(x):
    return numerix.maximum(x, 1e-8)

def safe_exp(x):
    return numerix.exp(numerix.clip(x, -100, 100))

def safe_div(numer, denom, eps=1e-12):
    return numer / numerix.maximum(denom, eps)

def radial_avg(Cvar, Nr, Nz):
    """Return radial average along reactor length (axis=0 is radial)."""
    return np.mean(Cvar.reshape(Nr, Nz), axis=0)
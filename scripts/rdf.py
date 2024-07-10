import os
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('ignore') # np log is annoying

class RDF:
    """
    Radial Distribution Function Class
    """
    def __init__(self,name,ID):
        self.name = name
        self.ID = ID
        self.r_data = None
        self.g_data = None
        self.G = None
        self.bounds = None

    def set_RDF(self,r,g):
        """
        Set r and g(r) of radial distributrion function
        """
        self.r_data = r
        self.g_data = g

    def get_r(self,g):
        """
        Given a g-value return r
        """
        return np.interp(g,self.r_data,self.g_data)

    def g(self,r):
        """
        Given a r-value return g
        """
        return np.interp(r,self.g_data,self.r_data)


    def get_G(self):
        """
        Galculate G(r)
        """
        kb = 8.31446261815324
        T = 300
        _sum = np.sum(self.g_data)
        self.G = -kb*T*np.log(self.g_data)/_sum


class Bounded_RDF(RDF):
    """
    Bounded Radial Distibution Class
    """
    def __init__(self,rdf):
        self.rdf = rdf
        self.r_data = None
        self.g_data = None
        self.bounds = [None,None]

    def set_bounds(self):
        """
        Set the bounds of the Bounded RDF
        """
        self.bounds[0] = self.find_min_r()
        self.bounds[1] = self.find_max_r(1.1)
        self.set_RDF()

    def find_min_r(self):
        """
        Find the smallest r values FROM DATA such that all g(r) values are non-zero after r
        """
        r_loc = np.where([self.rdf.g_data == 0])[1][-1]

        return r_loc

    def find_max_r(self,g):
        """
        Find the smallest r values FROM DATA such that all g(r) values are non-zero after r
        """
        find_r = g - self.rdf.g_data
        r_loc  = np.where([find_r < 0])[1][0]

        return r_loc

    def set_RDF(self):
        """
        Set the Bounds of the Radial Distribution Function
        """    
        self.r_data = self.rdf.r_data[self.bounds[0]:self.bounds[1]]
        self.g_data = self.rdf.g_data[self.bounds[0]:self.bounds[1]]

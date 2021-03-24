import numpy as np
import copy
import random

class Genome:
    def __init__(self):
        self.gene_size = 64
        self._genes = np.zeros(self.gene_size, dtype=bool)

    @staticmethod
    def set_mutation_rate(num):
        Genome._mutation_rate = num
    
    @staticmethod
    def set_bad_threshold(num):
	    Genome.bad_threshold = num
    
    def mutate(self):
	    indices = np.arange(0,self.gene_size)
	    np.random.shuffle(indices)
	
	    for k in indices[:Genome._mutation_rate]:
	        self._genes[k] = not self._genes[k]
	
    def count_bad(self, age):
        return sum(self._genes[:age])

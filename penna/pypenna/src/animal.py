from src.genome import Genome
import random
import copy

class Animal:
    def __init__(self, genome=Genome(), age = 0):
        self._genome = genome
        self._age = age
        self._pregnant = False
	
    @staticmethod
    def set_reproduction_age(num):
        Animal._reproduction_age = num
    
    @staticmethod
    def set_prob_to_get_pregnant(num):
        Animal._prob_to_get_pregnant = num
    
    def is_pregnant(self):
        return self._pregnant
    
    def age(self):
        return self._age
    
    def is_dead(self):
        return (self._genome.count_bad(self._age)>Genome.bad_threshold) \
                or (self._age > self._genome.gene_size)

    def give_birth(self):
        assert(self.is_pregnant())
        baby_genes = copy.deepcopy(self._genome)
        baby_genes.mutate()
        return Animal(baby_genes)
   
    def advance_by_one_year(self):
        self._age += 1
        if ( self._age > Animal._reproduction_age ) and \
            ( not self._pregnant ):
            if (Animal._prob_to_get_pregnant > random.random()):
                self._pregnant = True

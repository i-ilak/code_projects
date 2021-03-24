from src.animal import Animal
from src.genome import Genome
import random
import copy

class Population:
    def __init__(self, how_many, nmax=10000):
        self._pop = [Animal() for x in range(how_many)]
        self._nmax = nmax
   
    def size(self):
        return len(self._pop)

    def step(self):
        for animal in self._pop:
            animal.advance_by_one_year()

        self._pop[:] = [animal for animal in self._pop 
                        if not animal.is_dead()
                        if not random.random() < self.size()/self._nmax
                        ] 
        babies = [animal.give_birth() for animal in self._pop if animal.is_pregnant()]
        
        self._pop += babies
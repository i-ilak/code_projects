from src.population import Population
import random


def FishingPrediacte(animal, f_rate, f_age):
    return animal.age() > f_age and f_rate > random.random()


class FishingPopulation(Population):
    def __init__(self, how_many, nmax=10000):
        super().__init__(how_many, nmax)
        self._fishing_rate = 0
        self._fishing_age = 0

    def change_fishing(self, n_frate, n_fage):
        self._fishing_rate = n_frate
        self._fishing_age = n_fage

    def step(self):
        super().step()
        self._pop[:] = [animal for animal in self._pop 
                        if not FishingPrediacte(animal, 
                            self._fishing_rate, self._fishing_age)]
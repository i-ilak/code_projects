from src.fishing_pop import FishingPopulation
from src.animal import Animal
from src.genome import Genome
from tqdm import tqdm

if __name__ == "__main__":
    Genome.set_bad_threshold(2)
    Genome.set_mutation_rate(2)
    Animal.set_prob_to_get_pregnant(1)
    Animal.set_reproduction_age(8)
    
    fish = FishingPopulation(300)
    
    file = open("data.txt", "w")
    file.flush()
    for k in tqdm(range(1500)):
        fish.step()
        file.write("{}\t{}\n".format(k+1,fish.size()))
        if(k == 600):
            fish.change_fishing(0.17, 8)
        if(k == 1000):
            fish.change_fishing(0.5, 0)
    file.close()


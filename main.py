# IMPORTS
# import time
import random
import json
from CIEA import CIEA

### CONFIGURATION SPECS ###
# Load the JSON config file
configuration_Data = json.load(open("config_beckmann.json"))
# importing and parsing the configuration files
runs = configuration_Data["runs"]
offspringNumber = configuration_Data["offspringNumber"]
popSize = configuration_Data["popSize"]
generations = configuration_Data["generations"]
# 0 is perfect toxins, 1 is de-novo random, 2 is founder mutation instantiation
instantiationToggle = configuration_Data["instantiationToggle"]
# 0 is self evolvable and anything else is set.
mutationRate = configuration_Data["mutationRate"]
# 0 is off, anything else is number of generations with same avg fitness
terminationCondition = configuration_Data["terminationCondition"]
# setting the random seed for reproducibility
# if JSON is set to seed = 0 it will randomize the run
if configuration_Data["seed"] == 0:
    random.seed(random.randint(0, 100000))
else:
    # otherwise it will run on a predetermined set seed.
    random.seed(configuration_Data["seed"])

# instantiate the CIEA class
newCIEA = CIEA(instantiationToggle, mutationRate, terminationCondition)

### RUNNING A SIMULATION ###
# Run a classic TA simulation
# Here I am sending in Mew (population size)
newCIEA.runClassicTA(generations, runs, offspringNumber, popSize, mutationRate)


# clocking checker
"""
    print("testing time of survival selection")
    time0 = time.perf_counter()
    #method goes here
    time1 = time.perf_counter()
    timeItTakes = time1 - time0
    print("It took this long")
    print(timeItTakes)
    print("\n")
"""

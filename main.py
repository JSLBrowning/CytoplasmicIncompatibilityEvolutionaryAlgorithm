"""main driver"""
#IMPORTS
import time
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
instantiationToggle = configuration_Data["instantiationToggle"] #0 is perfect toxins, 1 is de-novo random, 2 is founder mutation instantiation
mutationRate = configuration_Data["mutationRate"] #0 is self evolvable and anything else is set.
terminationCondition = configuration_Data["terminationCondition"] #0 is off, anything else is number of generations with same avg fitness
#setting the random seed for reproducibility
if configuration_Data["seed"] == 0: #if JSON is set to seed = 0 it will randomize the run
    random.seed(random.randint(0, 100000))
else:
    random.seed(configuration_Data["seed"]) #otherwise it will run on a predetermined set seed.

#instantiate the CIEA class
newCIEA = CIEA(instantiationToggle, mutationRate, terminationCondition)

### RUNNING A SIMULATION ###
# Run a classic TA simulation
newCIEA.runClassicTA(generations, runs, offspringNumber, popSize, mutationRate) # Here I am sending in Mew (population size)















#clocking checker
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
"""This is an EA driver to drive the CI evolutionary simulation"""

# IMPORTS
from TA import TA
import random
import time
from Bio import pairwise2
from Bio.Seq import Seq


class CIEA:
    def __init__(self, instantiationToggle, mutationRate, terminationCondition):
        # CLASS VARIABLES
        self.mutationRate = mutationRate
        self.instantiationToggle = instantiationToggle
        self.TApopulation = []
        self.matingPool = []
        self.offspringPool = []
        self.avgFitness = -10000000000
        self.lastAvgFitness = -10000000000
        self.terminationCondition = terminationCondition
        self.terminationCounter = 0
        self.generationsCount = 0
        self.runi = 0

    # CLASS RESETS
    def resetEverything(self):
        self.resetmatingPool()
        self.resetOffspringPool()

    def resetTAPopulation(self):
        self.TApopulation = []
        self.avgFitness = -10000000000
        self.lastAvgFitness = -10000000000
        self.terminationCounter = 0
        self.generationsCount = 0

    def resetmatingPool(self):
        self.matingPool = []

    def resetOffspringPool(self):
        self.offspringPool = []

    # CLASS METHODS
    def generateSingleTA(self):
        newTA = TA(self.instantiationToggle, self.mutationRate)
        return newTA

    def generateTAPopulation(self, mew):  # mew is population size
        if self.instantiationToggle == 2:
            founder = self.generateSingleTA()
            self.TApopulation.append(founder)
            for i in range(1, mew):
                childTA = self.generateSingleTA()
                childTA.setToxinSchema(founder.toxinSchema)
                childTA.setAntidoteSchema(founder.antidoteSchema)
                self.mutateTA(childTA)
                self.TApopulation.append(childTA)
            for TA in self.TApopulation:
                TA.re_evaluateAllFitness()
            self.sortTAPopulation()
        else:
            for i in range(0, mew):
                self.TApopulation.append(self.generateSingleTA())
            self.sortTAPopulation()

    # Sorts the TA population based on fitness
    def sortSpecificTAPopulation(self, population):
        population.sort(reverse=True, key=lambda x: x.combinedTAFitness)

    # Sorts the TA population based on fitness
    def sortTAPopulation(self):
        self.TApopulation.sort(reverse=True, key=lambda x: x.combinedTAFitness)

    # ////METHODS FOR PARENT SELECTION
    # 0 1 toggle triggers sort by local fitness, or universal fitness
    # you need double the numberOffspring - of mates
    def PSTruncation(self, numberOffspring, population, toggle):
        if toggle == 0:
            self.sortSpecificTAPopulation(population)
        elif toggle == 1:
            self.sortUniversalTAPopulation(population)
        numberMates = numberOffspring*2
        for i in range(0, numberMates):
            self.matingPool.append(population[i])

    def PSKTournament(self, numberOffspring, k):
        numberMates = numberOffspring*2
        while len(self.matingPool) < numberMates:
            picks = random.sample(self.TApopulation, k=k)
            self.sortSpecificTAPopulation(picks)
            theOneToMate = picks[0]
            self.matingPool.append(theOneToMate)

    def PSFPS(self, numberOffspring, population):
        numberOfParentsNeeded = numberOffspring*2
        fitness_sum = sum([TA.combinedTAFitness for TA in population])
        w = [TA.combinedTAFitness / fitness_sum for TA in population]
        fpsPicks = random.choices(
            population, weights=w, k=numberOfParentsNeeded)
        for TA in fpsPicks:
            self.matingPool.append(TA)

    # A method for recombination using single-point crossover; input are two parent TA's; creates a child
    def variable_point_crossover_recombination(self, mom, dad):
        # Determine which toxin schema is smaller and which is larger.
        smallest_toxin, largest_toxin = dad.toxinSchema, mom.toxinSchema
        if len(mom.toxinSchema) < len(dad.toxinSchema):
            smallest_toxin, largest_toxin = mom.toxinSchema, dad.toxinSchema

        # Same, but for antidote schema.
        smallest_antidote, largest_antidote = dad.AntidoteSchema, mom.AntidoteSchema
        if len(mom.AntidoteSchema) < len(dad.AntidoteSchema):
            smallest_antidote, largest_antidote = mom.AntidoteSchema, dad.AntidoteSchema

        # Determine number of crossover points.
        crossover_points = random.randit(1, 6)

        # Randomly generate acceptable crossover points (then sort them).
        toxin_crossover_points = sorted(random.sample(
            range(1, len(smallest_toxin) - 1), crossover_points))
        antidote_crossover_points = sorted(random.sample(
            range(1, len(smallest_antidote) - 1), crossover_points))

        # Create a new toxin from the two parents.
        new_toxin = []
        for i in range(len(toxin_crossover_points)):
            if i % 2 == 0:
                new_toxin += smallest_toxin[toxin_crossover_points[i]:toxin_crossover_points[i + 1]]
            else:
                new_toxin += largest_toxin[toxin_crossover_points[i]:toxin_crossover_points[i + 1]]

        # Create a new antidote from the two parents.
        new_antidote = []
        for i in range(len(antidote_crossover_points)):
            if i % 2 == 0:
                new_antidote += smallest_antidote[antidote_crossover_points[i]:antidote_crossover_points[i + 1]]
            else:
                new_antidote += largest_antidote[antidote_crossover_points[i]:antidote_crossover_points[i + 1]]

        # Use new toxin and antidote to create a new TA.
        childTA = TA(self.instantiationToggle, self.mutationRate)
        childTA.setToxinSchema(new_toxin)
        childTA.setAntidoteSchema(new_antidote)

        # Optional mutation rate assignment.
        coin_flip = random.randint(0, 2)
        if coin_flip == 0:
            childTA.setMutationRate(dad.mutationRate)
        elif coin_flip == 1:
            childTA.setMutationRate(mom.mutationRate)
        elif coin_flip == 2 and self.mutationRate == 0:
            # Max mutation rate is 100.
            childTA.setMutationRate(random.randint(1, 100))

        self.offspringPool.append(childTA)

    # drives recombination (picking the most fit to mate with the next most fit, down the line promiscuously)

    def recombinationDriver(self, numberOffspring):
        i = 0
        while i < numberOffspring:
            coinFlip1 = random.randint(1, 1)
            if coinFlip1 == 0:
                # elitist sex (somewhat incestuous).
                self.variable_point_crossover_recombination(
                    self.matingPool[i], self.matingPool[i + 1])
                i += 1
            if coinFlip1 == 1:
                self.variable_point_crossover_recombination(self.matingPool[i], random.choice(
                    self.matingPool))  # this method is less incestious but still elitist
                i += 1
            else:
                self.variable_point_crossover_recombination(self.matingPool[i], random.choice(
                    self.TApopulation))  # promisciuos sex with a less fit individual
                i += 1

    # MUTATION OPERATOR FUNCTIONS FOR mutateTA().
    # Bit flip mutation.

    def bit_flip_mutation(self):
        if random.randint(0, 1) == 0:
            mutationIndex = random.randint(0, len(self.toxinSchema) - 1)
            newSchema = self.toxinSchema[0: mutationIndex] + random.choice(
                self.DNABases) + self.toxinSchema[mutationIndex + 1: len(self.toxinSchema)]
            if len(newSchema) <= self.maxToxinLength:
                self.setToxinSchema(newSchema)
        else:
            mutationIndex = random.randint(0, len(self.antidoteSchema) - 1)
            newSchema = self.antidoteSchema[0: mutationIndex] + random.choice(
                self.DNABases) + self.antidoteSchema[mutationIndex + 1: len(self.antidoteSchema)]
            if len(newSchema) <= self.maxAntidoteLength:
                self.setAntidoteSchema(newSchema)

    # Insertion mutation.

    def insertion_mutation(self):
        toxin_or_antidote = random.randint(0, 1)  # decides toxin or antidote
        # decides how much to insert (this is a magic number here) 12 means max 4 codons essentially
        insertion_length = random.randint(1, 12)
        if toxin_or_antidote == 0 and len(self.toxinSchema) != self.maxToxinLength:
            mutationIndex = random.randint(0, len(self.toxinSchema) - 1)
            sectionToAdd = ""
            for i in range(0, insertion_length):
                sectionToAdd = sectionToAdd + random.choice(self.DNABases)
            newSchema = self.toxinSchema[0: mutationIndex] + sectionToAdd + \
                self.toxinSchema[mutationIndex: len(self.toxinSchema)]
            if len(newSchema) <= self.maxAntidoteLength and len(newSchema) >= self.minAntidoteLength:
                self.setToxinSchema(newSchema)
        else:
            mutationIndex = random.randint(0, len(self.antidoteSchema) - 1)
            sectionToAdd = ""
            for i in range(0, insertion_length):
                sectionToAdd = sectionToAdd + random.choice(self.DNABases)
            newSchema = self.antidoteSchema[0: mutationIndex] + sectionToAdd + \
                self.antidoteSchema[mutationIndex: len(self.antidoteSchema)]
            if len(newSchema) <= self.maxAntidoteLength and len(newSchema) >= self.minAntidoteLength:
                self.setAntidoteSchema(newSchema)

    # Deletion mutation.

    def deletion_mutation(self):
        toxin_or_antidote = random.randint(0, 1)  # decides toxin or antidote
        deletion_length = random.randint(1, 12)  # decides how much to delete
        if toxin_or_antidote == 0:
            mutationIndex = random.randint(1, len(self.toxinSchema) - 1)
            if mutationIndex - deletion_length > 0:
                newSchema = self.toxinSchema[0: mutationIndex - deletion_length] + \
                    self.toxinSchema[mutationIndex: len(self.toxinSchema)]
                if len(newSchema) >= self.minToxinLength and len(newSchema) <= self.maxToxinLength:
                    self.setToxinSchema(newSchema)
        elif toxin_or_antidote == 1 and len(self.antidoteSchema) != self.minAntidoteLength:
            mutationIndex = random.randint(1, len(self.antidoteSchema) - 1)
            if mutationIndex - deletion_length > 0:
                newSchema = self.antidoteSchema[0: mutationIndex - deletion_length] + \
                    self.antidoteSchema[mutationIndex: len(
                        self.antidoteSchema)]
                if len(newSchema) >= self.minAntidoteLength:
                    self.setAntidoteSchema(newSchema)

    # Mutates one TA by a single bit flip, an insertion, or a deletion.

    def mutateTA(self, child):
        for _ in range(0, child.mutationRate):
            mutation_operator = random.randint(0, 3)  # 0 = do nothing
            if mutation_operator == 1:  # bit flip
                child.bit_flip_mutation()
            if mutation_operator == 2:  # insertion
                child.insertion_mutation()
            if mutation_operator == 3:  # deletion
                child.deletion_mutation()

    def mutationDriver(self, timesToMutate):
        for i in range(0, timesToMutate):
            for children in self.offspringPool:
                self.mutateTA(children)
        for children in self.offspringPool:
            children.re_evaluateAllFitness()

    # ////SURVIVOR SELECTION METHODS ////////

    def SSTruncation(self, immigrants):
        # accounts for immigration ... this would be buggy if you change stuff
        count = 0 + immigrants
        for children in self.offspringPool:
            self.TApopulation.append(children)
            count += 1
        self.sortTAPopulation()
        for int in range(0, count):
            self.TApopulation.pop()

    def SSKTourny(self, mew, k):
        # select k, 1 wins, decide what to do
        for children in self.offspringPool:
            self.TApopulation.append(children)
        if k < mew:
            newPopulation = []
            while len(self.TApopulation) > k - 1 and len(newPopulation) < mew:
                picks = random.sample(self.TApopulation, k=k)
                self.sortSpecificTAPopulation(picks)
                theOneToSave = picks[0]
                self.TApopulation.remove(theOneToSave)
                newPopulation.append(theOneToSave)
        self.sortSpecificTAPopulation(newPopulation)
        self.TApopulation = newPopulation

    # READOUT METHODS
    def calculateAvgFitness(self):
        fitnessSum = sum([TA.combinedTAFitness for TA in self.TApopulation])
        avg = fitnessSum / len(self.TApopulation)
        return avg

    def calculateAvgBindingFitness(self):
        fitnessSum = sum([TA.localBindingFitness for TA in self.TApopulation])
        avg = fitnessSum / len(self.TApopulation)
        return avg

    def calculateAvgDUBFitness(self):
        fitnessSum = sum(
            [TA.bestCatalyticDUBScore for TA in self.TApopulation])
        avg = fitnessSum / len(self.TApopulation)
        return avg

    def calculateAvgNucFitness(self):
        fitnessSum = sum(
            [TA.bestCatalyticNucScore for TA in self.TApopulation])
        avg = fitnessSum / len(self.TApopulation)
        return avg

    def calculateAvgToxinLength(self):
        toxinLength = sum([len(TA.toxinSchema) for TA in self.TApopulation])
        avg = toxinLength / len(self.TApopulation)
        return avg

    def calculateAvgAntidoteLength(self):
        antidoteLength = sum([len(TA.antidoteSchema)
                             for TA in self.TApopulation])
        avg = antidoteLength / len(self.TApopulation)
        return avg

    def calculateAvgToxinTranslationLength(self):
        toxinLength = sum([len(TA.toxinTranslation)
                          for TA in self.TApopulation])
        avg = toxinLength / len(self.TApopulation)
        return avg

    def calculateAvgAntidoteTranslationLength(self):
        antidoteLength = sum([len(TA.antidoteTranslation)
                             for TA in self.TApopulation])
        avg = antidoteLength / len(self.TApopulation)
        return avg

    def calculateAvgTAMutationRate(self):
        toxinLength = sum([TA.mutationRate for TA in self.TApopulation])
        avg = toxinLength / len(self.TApopulation)
        return avg

    def calculateAvg_NLS_Fitness(self):
        fitnessSum = sum([TA.bestNLSScore for TA in self.TApopulation])
        avg = fitnessSum / len(self.TApopulation)
        return avg

    def calculateAvg_T4SS_Fitness(self):
        fitnessSum = sum([TA.bestTypeIVSSScore for TA in self.TApopulation])
        avg = fitnessSum / len(self.TApopulation)
        return avg

    def calculateAvgNLSLocation(self):
        SiteSum = sum([TA.NLSSiteLocation for TA in self.TApopulation])
        avg = SiteSum / len(self.TApopulation)
        return avg

    def calculateAvgTypeIVLocation(self):
        SiteSum = sum([TA.TypeIVSiteLocation for TA in self.TApopulation])
        avg = SiteSum / len(self.TApopulation)
        return avg

    def calculateDiversityIndex(self):
        # this takes a sample of ten 10's and aligns them to each other, summing their scores (high score means low diversity)
        # pick a random sample of 10 individuals
        sample = random.choices(self.TApopulation, k=11)
        diversityIndex = 0
        for i in range(0, len(sample)-1):
            diversityIndex = diversityIndex + float(pairwise2.align.globalxs(
                sample[i].toxinTranslation, sample[i+1].toxinTranslation, -1, -.1, penalize_end_gaps=False, score_only=True))
            i += 1
        return diversityIndex/10  # this is now normalizing with length in the denominator

    def immigrationTA(self, emigrationPopSize):
        self.generateTAPopulation(emigrationPopSize)

    # METHODS for SIMULATION
    def runClassicTA(self, generations, runs, offspringNumber, popSize, mutationRate):
        print("initiating Classic TA simulation...\n")
        self.reset_Results_Log()
        self.reset_Solutions_Log()

        runCount = 0
        while runCount < runs:
            self.generateTAPopulation(popSize)
            while self.generationsCount < generations:
                start = time.time()
                self.PSKTournament(offspringNumber, 5)  # ktourny number
                # self.PSFPS(offspringNumber, self.TApopulation)
                # self.PSTruncation(offspringNumber, self.TApopulation, 0)
                self.recombinationDriver(offspringNumber)
                self.mutationDriver(mutationRate)
                if self.instantiationToggle != 2:
                    # this method is immigrating into the population some new blood
                    self.immigrationTA(popSize//100)
                # self.SSKTourny(popSize, 5)
                if self.instantiationToggle != 2:
                    # the parameter here is accounting for immigration
                    self.SSTruncation(popSize//100)
                else:
                    # the parameter here is accounting for immigration
                    self.SSTruncation(0)
                self.resetEverything()
                print("\nbest TA schema: \n")
                self.TApopulation[0].toString()
                print("population size" + str(len(self.TApopulation)))
                print("Top 20 Toxins: \n")
                for i in range(0, 20):
                    print(self.TApopulation[i].toxinTranslation)

                # things to measure:
                highestTAFitness_HTF = " HTF_" + \
                    str(self.TApopulation[0].combinedTAFitness)
                avgBindingFitness_ABF = " ABF_" + \
                    str(self.calculateAvgBindingFitness())
                avgDUBFitness_ADF = " ADF_" + \
                    str(self.calculateAvgDUBFitness())
                avgNucFitness_ANF = " ANF_" + \
                    str(self.calculateAvgNucFitness())
                avgTAfitness_ATF = " ATF_" + str(self.calculateAvgFitness())
                avgToxinLength_ATL = " ATL_" + \
                    str(self.calculateAvgToxinLength())
                avgToxinAALength_ATAL = " ATAL_" + \
                    str(self.calculateAvgToxinTranslationLength())
                avgAntidoteLength_AAL = " AAL_" + \
                    str(self.calculateAvgAntidoteLength())
                avgAntidoteAALength_AAAL = " AAAL_" + \
                    str(self.calculateAvgAntidoteTranslationLength())
                avgTAMutationRate_ATMR = " ATMR_" + \
                    str(self.calculateAvgTAMutationRate())
                avgNLSSITELocation = " NLSL_" + \
                    str(self.calculateAvgNLSLocation())
                avgTypeIVSITELocation = " TYPL_" + \
                    str(self.calculateAvgTypeIVLocation())
                avgNLSFitness_ANLSF = " ANLSF_" + \
                    str(self.calculateAvg_NLS_Fitness())
                avgT4SSFitness_AT4F = " AT4F_" + \
                    str(self.calculateAvg_T4SS_Fitness())
                diversityIndex = " DI_" + \
                    str(1-(self.calculateDiversityIndex() /
                        self.calculateAvgToxinTranslationLength()))
                stringResults = highestTAFitness_HTF + avgBindingFitness_ABF + avgDUBFitness_ADF + avgNucFitness_ANF + avgTAfitness_ATF + avgToxinLength_ATL + avgToxinAALength_ATAL + \
                    avgAntidoteLength_AAL + avgAntidoteAALength_AAAL + avgTAMutationRate_ATMR + avgNLSSITELocation + \
                    avgTypeIVSITELocation + avgNLSFitness_ANLSF + \
                    avgT4SSFitness_AT4F + diversityIndex
                self.write_Results_Log(stringResults)

                print("highest TA fitness: " + highestTAFitness_HTF)
                print("average Toxin AA Length: " + avgToxinAALength_ATAL)
                print("average antidote AA Length: " + avgAntidoteAALength_AAAL)
                print("avg binding fitness: " + avgBindingFitness_ABF)
                print("avg DUB fitness: " + avgDUBFitness_ADF)
                print("avg Nuc fitness: " + avgNucFitness_ANF)
                print("avg TA fitness: " + avgTAfitness_ATF)
                print("avg TA mutation rate: " + avgTAMutationRate_ATMR)
                print("avg NLS Site: " + avgNLSSITELocation)
                print("avg Type IV Site: " + avgTypeIVSITELocation)
                print("avg NLS Fitness: " + avgNLSFitness_ANLSF)
                print("avg Type IV Fitness: " + avgT4SSFitness_AT4F)
                print("diversity index: " + diversityIndex)
                self.generationsCount += 1
                end = time.time()
                print("this generation took this many seconds: " + str(end - start))
                print("Generation: " + str(self.generationsCount))
                print("termination counter = " + str(self.terminationCounter))
                print("pop number" + str(len(self.TApopulation)))
                if self.checkTerminationConditions() == 0:  # terminate if fitness has not changed for # generations
                    break
            self.write_Solutions_Log(self.TApopulation[0].toString(
            ) + "\n simulation ended after: " + str(self.generationsCount) + " generations")
            self.resetTAPopulation()
            runCount += 1
            self.runi = runCount
            self.separate_Results_Logs()

    # ////LOGGING METHODS////////
    # a function to clear the text file for new results
    def reset_Results_Log(self):
        report_Results_log = open("resultsLog.txt", "w")
        report_Results_log.write("Run 0\n")
        report_Results_log.close()

    # a function to deliminate run data
    def separate_Results_Logs(self):
        print("separate results log")
        report_Results_log = open("resultsLog.txt", "a")
        report_Results_log.write("\nRun " + str(self.runi) + "\n")
        report_Results_log.close()

    # a function to reset the results log
    def reset_Solutions_Log(self):
        report_Results_log = open("SolutionLog.txt", "w")
        report_Results_log.write("Run 0\n")
        report_Results_log.close()

    # a function to separate logs of the best solutions in the solutions file
    def separate_Solution_Logs(self):
        report_Solutions_log = open("SolutionLog.txt", "a")
        report_Solutions_log.write("\nRun " + str(self.runi) + "\n")
        report_Solutions_log.close()

    # logging the results information function.
    def write_Results_Log(self, stringResult):
        report_Results_log = open("resultsLog.txt", "a")
        report_Results_log.write(stringResult + "\n")
        report_Results_log.close()

    # logging the solutions information function.
    def write_Solutions_Log(self, string):
        report_Solutions_log = open("SolutionLog.txt", "a")
        report_Solutions_log.write(string + "\n")
        report_Solutions_log.close()

    def checkTerminationConditions(self):
        if self.terminationCondition == 0:
            return 1
        else:
            self.lastAvgFitness = self.avgFitness
            self.avgFitness = self.calculateAvgFitness()
            if abs(self.avgFitness - self.lastAvgFitness) <= 0.001:
                self.terminationCounter += 1
                if self.terminationCounter > self.terminationCondition:  # fitness has not changed for 25 generations
                    return 0
                else:
                    return 1
            else:
                self.terminationCounter = 0
                return 1

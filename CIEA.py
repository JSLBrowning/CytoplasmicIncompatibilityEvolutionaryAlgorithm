"""This is an EA driver to drive the CI evolutionary simulation"""

#IMPORTS
from TA import TA
import random
import time
from Bio import pairwise2
from Bio.Seq import Seq

class CIEA:
    def __init__(self, instantiationToggle, mutationRate, terminationCondition):
        #CLASS VARIABLES
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

    #CLASS RESETS
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

    #CLASS METHODS
    def generateSingleTA(self):
        newTA = TA(self.instantiationToggle, self.mutationRate)
        return newTA

    def generateTAPopulation(self, mew): #mew is population size
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

    #Sorts the TA population based on fitness
    def sortSpecificTAPopulation(self, population):
        population.sort(reverse=True, key=lambda x: x.combinedTAFitness)

    #Sorts the TA population based on fitness
    def sortTAPopulation(self):
        self.TApopulation.sort(reverse=True, key=lambda x: x.combinedTAFitness)

    #////METHODS FOR PARENT SELECTION
    #0 1 toggle triggers sort by local fitness, or universal fitness
    def PSTruncation(self, numberOffspring, population, toggle): # you need double the numberOffspring - of mates
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
        fpsPicks = random.choices(population, weights=w, k=numberOfParentsNeeded)
        for TA in fpsPicks:
            self.matingPool.append(TA)

    #A method for recombination using single-point crossover; input are two parent TA's; creates a child
    def recombineTA(self, mom, dad):
        #figures out which is smaller and larger schema
        if len(mom.toxinSchema) < len(dad.toxinSchema):
            smallestToxin = mom.toxinSchema
            largestToxin = dad.toxinSchema
        else:
            smallestToxin = dad.toxinSchema
            largestToxin = mom.toxinSchema
        if len(mom.antidoteSchema) < len(dad.antidoteSchema):
            smallestAntidote = mom.antidoteSchema
            largestAntidote = dad.antidoteSchema
        else:
            smallestAntidote = dad.antidoteSchema
            largestAntidote = mom.antidoteSchema

        coinFlip = random.randint(3,3) #right now its only doing double point crossover
        if coinFlip == 0: #single point crossover
            toxinCrossoverPoint = random.randint(0, len(smallestToxin) - 1)
            antidoteCrossoverPoint = random.randint(0, len(smallestAntidote) - 1)
            newToxin = smallestToxin[0: toxinCrossoverPoint] + largestToxin[toxinCrossoverPoint: len(largestToxin)]
            newAntidote = smallestAntidote[0: antidoteCrossoverPoint] + largestAntidote[antidoteCrossoverPoint: len(largestAntidote)]
            childTA = TA(self.instantiationToggle, self.mutationRate)
            childTA.setToxinSchema(newToxin)
            childTA.setAntidoteSchema(newAntidote)
        elif coinFlip == 1: #double point crossover
            toxinCrossoverPoint = random.randint(0, len(smallestToxin) - 1)
            toxinCrossoverPoint2 = random.randint(toxinCrossoverPoint, len(smallestToxin) - 1)
            antidoteCrossoverPoint = random.randint(0, len(smallestAntidote) - 1)
            antidoteCrossoverPoint2 = random.randint(antidoteCrossoverPoint, len(smallestAntidote) - 1)

            newToxin = smallestToxin[0: toxinCrossoverPoint] + largestToxin[toxinCrossoverPoint: toxinCrossoverPoint2] + smallestToxin[toxinCrossoverPoint2: len(largestToxin)]
            newAntidote = smallestAntidote[0: antidoteCrossoverPoint] + largestAntidote[antidoteCrossoverPoint: antidoteCrossoverPoint2] + smallestAntidote[antidoteCrossoverPoint2: len(largestAntidote)]
            childTA = TA(self.instantiationToggle, self.mutationRate)
            childTA.setToxinSchema(newToxin)
            childTA.setAntidoteSchema(newAntidote)
        elif coinFlip == 2: #triple point crossover
            toxinCrossoverPoint = random.randint(0, len(smallestToxin) - 1)
            toxinCrossoverPoint2 = random.randint(toxinCrossoverPoint, len(smallestToxin) - 1)
            toxinCrossoverPoint3 = random.randint(toxinCrossoverPoint2, len(smallestToxin) - 1)
            antidoteCrossoverPoint = random.randint(0, len(smallestAntidote) - 1)
            antidoteCrossoverPoint2 = random.randint(antidoteCrossoverPoint, len(smallestAntidote) - 1)
            antidoteCrossoverPoint3 = random.randint(antidoteCrossoverPoint2, len(smallestAntidote) - 1)

            newToxin = smallestToxin[0: toxinCrossoverPoint] + largestToxin[toxinCrossoverPoint: toxinCrossoverPoint2] + smallestToxin[toxinCrossoverPoint2: toxinCrossoverPoint3] + largestToxin[toxinCrossoverPoint3: len(largestToxin)]
            newAntidote = smallestAntidote[0: antidoteCrossoverPoint] + largestAntidote[antidoteCrossoverPoint: antidoteCrossoverPoint2] + smallestAntidote[antidoteCrossoverPoint2: antidoteCrossoverPoint3]+ largestAntidote[antidoteCrossoverPoint3: len(largestAntidote)]
            childTA = TA(self.instantiationToggle, self.mutationRate)
            childTA.setToxinSchema(newToxin)
            childTA.setAntidoteSchema(newAntidote)
        elif coinFlip == 3: #4 point crossover
            toxinCrossoverPoint = random.randint(0, len(smallestToxin) - 1)
            toxinCrossoverPoint2 = random.randint(toxinCrossoverPoint, len(smallestToxin) - 1)
            toxinCrossoverPoint3 = random.randint(toxinCrossoverPoint2, len(smallestToxin) - 1)
            toxinCrossoverPoint4 = random.randint(toxinCrossoverPoint3, len(smallestToxin) - 1)
            antidoteCrossoverPoint = random.randint(0, len(smallestAntidote) - 1)
            antidoteCrossoverPoint2 = random.randint(antidoteCrossoverPoint, len(smallestAntidote) - 1)
            antidoteCrossoverPoint3 = random.randint(antidoteCrossoverPoint2, len(smallestAntidote) - 1)
            antidoteCrossoverPoint4 = random.randint(antidoteCrossoverPoint2, len(smallestAntidote) - 1)

            newToxin = smallestToxin[0: toxinCrossoverPoint] + largestToxin[toxinCrossoverPoint: toxinCrossoverPoint2] + smallestToxin[toxinCrossoverPoint2: toxinCrossoverPoint3] + largestToxin[toxinCrossoverPoint3: toxinCrossoverPoint4] + smallestToxin[toxinCrossoverPoint4: len(largestToxin)]
            newAntidote = smallestAntidote[0: antidoteCrossoverPoint] + largestAntidote[antidoteCrossoverPoint: antidoteCrossoverPoint2] + smallestAntidote[antidoteCrossoverPoint2: antidoteCrossoverPoint3]+ largestAntidote[antidoteCrossoverPoint3: antidoteCrossoverPoint4] + smallestAntidote[antidoteCrossoverPoint4: len(largestAntidote)]
            childTA = TA(self.instantiationToggle, self.mutationRate)
            childTA.setToxinSchema(newToxin)
            childTA.setAntidoteSchema(newAntidote)
        elif coinFlip == 4: #6 point crossover
            toxinCrossoverPoint = random.randint(0, len(smallestToxin) - 1)
            toxinCrossoverPoint2 = random.randint(toxinCrossoverPoint, len(smallestToxin) - 1)
            toxinCrossoverPoint3 = random.randint(toxinCrossoverPoint2, len(smallestToxin) - 1)
            toxinCrossoverPoint4 = random.randint(toxinCrossoverPoint3, len(smallestToxin) - 1)
            toxinCrossoverPoint5 = random.randint(toxinCrossoverPoint4, len(smallestToxin) - 1)
            toxinCrossoverPoint6 = random.randint(toxinCrossoverPoint5, len(smallestToxin) - 1)
            antidoteCrossoverPoint = random.randint(0, len(smallestAntidote) - 1)
            antidoteCrossoverPoint2 = random.randint(antidoteCrossoverPoint, len(smallestAntidote) - 1)
            antidoteCrossoverPoint3 = random.randint(antidoteCrossoverPoint2, len(smallestAntidote) - 1)
            antidoteCrossoverPoint4 = random.randint(antidoteCrossoverPoint2, len(smallestAntidote) - 1)
            antidoteCrossoverPoint5 = random.randint(antidoteCrossoverPoint2, len(smallestAntidote) - 1)
            antidoteCrossoverPoint6 = random.randint(antidoteCrossoverPoint2, len(smallestAntidote) - 1)

            newToxin = smallestToxin[0: toxinCrossoverPoint] + largestToxin[toxinCrossoverPoint: toxinCrossoverPoint2] + smallestToxin[toxinCrossoverPoint2: toxinCrossoverPoint3] + largestToxin[toxinCrossoverPoint3: toxinCrossoverPoint4] + smallestToxin[toxinCrossoverPoint4: toxinCrossoverPoint5] + largestToxin[toxinCrossoverPoint5: toxinCrossoverPoint6] + smallestToxin[toxinCrossoverPoint6: len(largestToxin)]
            newAntidote = smallestAntidote[0: antidoteCrossoverPoint] + largestAntidote[antidoteCrossoverPoint: antidoteCrossoverPoint2] + smallestAntidote[antidoteCrossoverPoint2: antidoteCrossoverPoint3] + largestAntidote[antidoteCrossoverPoint3: antidoteCrossoverPoint4] + smallestAntidote[antidoteCrossoverPoint4: antidoteCrossoverPoint5] + largestAntidote[antidoteCrossoverPoint5: antidoteCrossoverPoint6] + smallestAntidote[antidoteCrossoverPoint6: len(largestAntidote)]
            childTA = TA(self.instantiationToggle, self.mutationRate)
            childTA.setToxinSchema(newToxin)
            childTA.setAntidoteSchema(newAntidote)

        #optional mutation rate section
        coinFlip = random.randint(0,2)
        if coinFlip == 0:
            childTA.setMutationRate(dad.mutationRate)
        elif coinFlip == 1:
            childTA.setMutationRate(mom.mutationRate)
        elif coinFlip == 2 and self.mutationRate == 0:
            childTA.setMutationRate(random.randint(1,100)) #max mutation rate
        #end optional section

        self.offspringPool.append(childTA)

    #drives recombination (picking the most fit to mate with the next most fit, down the line promiscuously)
    def recombinationDriver(self, numberOffspring):
        i = 0
        while i < numberOffspring:
            coinFlip1 = random.randint(1, 1)
            if coinFlip1 == 0:
                self.recombineTA(self.matingPool[i], self.matingPool[i + 1]) #elitist sex (somewhat incestuous).
                i+=1
            if coinFlip1 == 1:
                self.recombineTA(self.matingPool[i], random.choice(self.matingPool)) #this method is less incestious but still elitist
                i+=1
            else:
                self.recombineTA(self.matingPool[i], random.choice(self.TApopulation)) #promisciuos sex with a less fit individual
                i+=1

    #Mutates TA by a single bit flip in 1 or the other; with a chance for deletion or addition as well.
    def mutateTA(self, child):
        for i in range(0, child.mutationRate):
            coinFlip1 = random.randint(0, 3) #0 do nothing, 1 = bit flip, 2 = insert bit, 3 = deletion
            if coinFlip1 == 1: #bit flip
                coinFlip2 = random.randint(0, 1)
                if coinFlip2 == 0:
                    mutationIndex = random.randint(0, len(child.toxinSchema) - 1)
                    newSchema = child.toxinSchema[0: mutationIndex] + random.choice(child.DNABases) + child.toxinSchema[mutationIndex + 1: len(child.toxinSchema)]
                    if len(newSchema) <= child.maxToxinLength:
                        child.setToxinSchema(newSchema)
                else:
                    mutationIndex = random.randint(0, len(child.antidoteSchema) - 1)
                    newSchema = child.antidoteSchema[0: mutationIndex] + random.choice(child.DNABases) + child.antidoteSchema[mutationIndex + 1: len(child.antidoteSchema)]
                    if len(newSchema) <= child.maxAntidoteLength:
                        child.setAntidoteSchema(newSchema)

            if coinFlip1 == 2: # insertion
                coinFlip2 = random.randint(0, 1) #decides toxin or antidote
                coinFlip3 = random.randint(1, 12)  # decides how much to insert (this is a magic number here) 12 means max 4 codons essentially
                if coinFlip2 == 0 and len(child.toxinSchema) != child.maxToxinLength:
                    mutationIndex = random.randint(0, len(child.toxinSchema) - 1)
                    sectionToAdd = ""
                    for i in range(0, coinFlip3):
                        sectionToAdd = sectionToAdd + random.choice(child.DNABases)
                    newSchema = child.toxinSchema[0: mutationIndex] + sectionToAdd + child.toxinSchema[mutationIndex: len(child.toxinSchema)]
                    if len(newSchema) <= child.maxAntidoteLength and len(newSchema) >= child.minAntidoteLength:
                        child.setToxinSchema(newSchema)
                else:
                    mutationIndex = random.randint(0, len(child.antidoteSchema) - 1)
                    sectionToAdd = ""
                    for i in range(0, coinFlip3):
                        sectionToAdd = sectionToAdd + random.choice(child.DNABases)
                    newSchema = child.antidoteSchema[0: mutationIndex] + sectionToAdd + child.antidoteSchema[mutationIndex: len(child.antidoteSchema)]
                    if len(newSchema) <= child.maxAntidoteLength and len(newSchema) >= child.minAntidoteLength:
                        child.setAntidoteSchema(newSchema)

            if coinFlip1 == 3: # deletion
                coinFlip2 = random.randint(0, 1) #decides toxin or antidote
                coinFlip3 = random.randint(1, 12) #decides how much to delete
                if coinFlip2 == 0:
                    mutationIndex = random.randint(1, len(child.toxinSchema) - 1)
                    if mutationIndex - coinFlip3 > 0:
                        newSchema = child.toxinSchema[0: mutationIndex - coinFlip3] + child.toxinSchema[mutationIndex: len(child.toxinSchema)]
                        if len(newSchema) >= child.minToxinLength and len(newSchema) <= child.maxToxinLength:
                            child.setToxinSchema(newSchema)
                elif coinFlip2 == 1 and len(child.antidoteSchema) != child.minAntidoteLength:
                    mutationIndex = random.randint(1, len(child.antidoteSchema) - 1)
                    if mutationIndex - coinFlip3 > 0:
                        newSchema = child.antidoteSchema[0: mutationIndex - coinFlip3] + child.antidoteSchema[mutationIndex: len(child.antidoteSchema)]
                        if len(newSchema) >= child.minAntidoteLength:
                            child.setAntidoteSchema(newSchema)

    def mutationDriver(self, timesToMutate):
        for i in range(0, timesToMutate):
            for children in self.offspringPool:
                self.mutateTA(children)
        for children in self.offspringPool:
            children.re_evaluateAllFitness()


    #////SURVIVOR SELECTION METHODS ////////
    def SSTruncation(self, immigrants):
        count = 0 + immigrants #accounts for immigration ... this would be buggy if you change stuff
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

    #READOUT METHODS
    def calculateAvgFitness(self):
        fitnessSum = sum([TA.combinedTAFitness for TA in self.TApopulation])
        avg = fitnessSum / len(self.TApopulation)
        return avg

    def calculateAvgBindingFitness(self):
        fitnessSum = sum([TA.localBindingFitness for TA in self.TApopulation])
        avg = fitnessSum / len(self.TApopulation)
        return avg

    def calculateAvgDUBFitness(self):
        fitnessSum = sum([TA.bestCatalyticDUBScore for TA in self.TApopulation])
        avg = fitnessSum / len(self.TApopulation)
        return avg

    def calculateAvgNucFitness(self):
        fitnessSum = sum([TA.bestCatalyticNucScore for TA in self.TApopulation])
        avg = fitnessSum / len(self.TApopulation)
        return avg

    def calculateAvgToxinLength(self):
        toxinLength = sum([len(TA.toxinSchema) for TA in self.TApopulation])
        avg = toxinLength / len(self.TApopulation)
        return avg

    def calculateAvgAntidoteLength(self):
        antidoteLength = sum([len(TA.antidoteSchema) for TA in self.TApopulation])
        avg = antidoteLength / len(self.TApopulation)
        return avg

    def calculateAvgToxinTranslationLength(self):
        toxinLength = sum([len(TA.toxinTranslation) for TA in self.TApopulation])
        avg = toxinLength / len(self.TApopulation)
        return avg

    def calculateAvgAntidoteTranslationLength(self):
        antidoteLength = sum([len(TA.antidoteTranslation) for TA in self.TApopulation])
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
        #this takes a sample of ten 10's and aligns them to each other, summing their scores (high score means low diversity)
        sample = random.choices(self.TApopulation, k=11) # pick a random sample of 10 individuals
        diversityIndex = 0
        for i in range(0,len(sample)-1):
            diversityIndex = diversityIndex + float(pairwise2.align.globalxs(sample[i].toxinTranslation, sample[i+1].toxinTranslation, -1, -.1, penalize_end_gaps=False, score_only=True))
            i+=1
        return diversityIndex/10 # this is now normalizing with length in the denominator

    def immigrationTA(self, emigrationPopSize):
        self.generateTAPopulation(emigrationPopSize)

    #METHODS for SIMULATION
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
                #self.PSFPS(offspringNumber, self.TApopulation)
                #self.PSTruncation(offspringNumber, self.TApopulation, 0)
                self.recombinationDriver(offspringNumber)
                self.mutationDriver(mutationRate)
                if self.instantiationToggle != 2:
                    self.immigrationTA(popSize//100) #this method is immigrating into the population some new blood
                #self.SSKTourny(popSize, 5)
                if self.instantiationToggle != 2:
                    self.SSTruncation(popSize//100) #the parameter here is accounting for immigration
                else:
                    self.SSTruncation(0) #the parameter here is accounting for immigration
                self.resetEverything()
                print("\nbest TA schema: \n")
                self.TApopulation[0].toString()
                print("population size" + str(len(self.TApopulation)))
                print("Top 20 Toxins: \n")
                for i in range (0,20):
                    print(self.TApopulation[i].toxinTranslation)

                #things to measure:
                highestTAFitness_HTF = " HTF_" + str(self.TApopulation[0].combinedTAFitness)
                avgBindingFitness_ABF = " ABF_" + str(self.calculateAvgBindingFitness())
                avgDUBFitness_ADF = " ADF_" + str(self.calculateAvgDUBFitness())
                avgNucFitness_ANF = " ANF_" + str(self.calculateAvgNucFitness())
                avgTAfitness_ATF = " ATF_" + str(self.calculateAvgFitness())
                avgToxinLength_ATL = " ATL_" + str(self.calculateAvgToxinLength())
                avgToxinAALength_ATAL = " ATAL_" + str(self.calculateAvgToxinTranslationLength())
                avgAntidoteLength_AAL = " AAL_" + str(self.calculateAvgAntidoteLength())
                avgAntidoteAALength_AAAL = " AAAL_" + str(self.calculateAvgAntidoteTranslationLength())
                avgTAMutationRate_ATMR = " ATMR_" + str(self.calculateAvgTAMutationRate())
                avgNLSSITELocation = " NLSL_" + str(self.calculateAvgNLSLocation())
                avgTypeIVSITELocation = " TYPL_" + str(self.calculateAvgTypeIVLocation())
                avgNLSFitness_ANLSF = " ANLSF_" + str(self.calculateAvg_NLS_Fitness())
                avgT4SSFitness_AT4F = " AT4F_" + str(self.calculateAvg_T4SS_Fitness())
                diversityIndex = " DI_" + str(1-(self.calculateDiversityIndex()/self.calculateAvgToxinTranslationLength()))
                stringResults = highestTAFitness_HTF + avgBindingFitness_ABF + avgDUBFitness_ADF + avgNucFitness_ANF + avgTAfitness_ATF + avgToxinLength_ATL + avgToxinAALength_ATAL + avgAntidoteLength_AAL + avgAntidoteAALength_AAAL + avgTAMutationRate_ATMR + avgNLSSITELocation + avgTypeIVSITELocation + avgNLSFitness_ANLSF + avgT4SSFitness_AT4F + diversityIndex
                self.write_Results_Log(stringResults)

                print("highest TA fitness: " + highestTAFitness_HTF)
                print("average Toxin AA Length: " + avgToxinAALength_ATAL)
                print("average antidote AA Length: " + avgAntidoteAALength_AAAL)
                print("avg binding fitness: " + avgBindingFitness_ABF)
                print("avg DUB fitness: " + avgDUBFitness_ADF)
                print("avg Nuc fitness: " + avgNucFitness_ANF)
                print("avg TA fitness: " + avgTAfitness_ATF)
                print("avg TA mutation rate: " +  avgTAMutationRate_ATMR)
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
                if self.checkTerminationConditions() == 0: #terminate if fitness has not changed for # generations
                    break
            self.write_Solutions_Log(self.TApopulation[0].toString() + "\n simulation ended after: " + str(self.generationsCount) + " generations")
            self.resetTAPopulation()
            runCount += 1
            self.runi = runCount
            self.separate_Results_Logs()

    #////LOGGING METHODS////////
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
                if self.terminationCounter > self.terminationCondition: #fitness has not changed for 25 generations
                    return 0
                else:
                    return 1
            else:
                self.terminationCounter = 0
                return 1




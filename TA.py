"""
This class is a Toxin-Antidote (TA) module analog.
version 1.0 - Initiated 5.19.2021.
version 1.1 - Finalized Draft Framework 6.17.2021.
version 1.2 - Updated 2.11.2022.
version 1.3 - Updated and cleaned up on 3.8.2022.
Author: Dr., Professor, and Auburn Senator John Frederick Beckmann
"""

# IMPORTS
import random
from Bio import pairwise2
from Bio.Seq import Seq


class TA:
    def __init__(self, instantiationToggle, mutationRate):
        # CLASS VARIABLES
        self.instantiationToggle = instantiationToggle
        self.DNABases = ["A", "T", "G", "C"]
        self.codonsNoSTOP = ["TTT", "TTC", "TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "ATT", "ATC", "ATA", "ATG", "GTT", "GTC", "GTA", "GTG", "TCT", "TCC", "TCA", "TCG", "CCT", "CCC", "CCA", "CCG", "ACT", "ACC", "ACA", "ACG", "GCT",
                             "GCC", "GCA", "GCG", "TAT", "TAC", "CAT", "CAC", "CAA", "CAG", "AAT", "AAC", "AAA", "AAG", "GAT", "GAC", "GAA", "GAG", "TGT", "TGC", "TGG", "CGT", "CGC", "CGA", "CGG", "AGT", "AGC", "AGA", "AGG", "GGT", "GGC", "GGA", "GGG"]
        self.codons = ["TTT", "TTC", "TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "ATT", "ATC", "ATA", "ATG", "GTT", "GTC", "GTA", "GTG", "TCT", "TCC", "TCA", "TCG", "CCT", "CCC", "CCA", "CCG", "ACT", "ACC", "ACA", "ACG", "GCT", "GCC", "GCA", "GCG", "TAT", "TAC",
                       "CAT", "CAC", "CAA", "CAG", "AAT", "AAC", "AAA", "AAG", "GAT", "GAC", "GAA", "GAG", "TGT", "TGC", "TGG", "CGT", "CGC", "CGA", "CGG", "AGT", "AGC", "AGA", "AGG", "GGT", "GGC", "GGA", "GGG", "TAG", "TAA", "TGA"]  # "TAG", "TAA", "TGA" stop codons
        # seed words are for instantiating with sub toxin words to make the evolution faster.
        self.seedWords = ["CCGATTATTATTGAACTGAA", "GATCTGGTGCTG", "CCGATTGGCCTGGAACTGAAA",
                          "CATTGGGTGACCCTGGTGATT", "TATTATGCGGATAGCCTG", "CAGCAGGCGGATGGCGCGGCGTGCGGC"]
        self.positiveAA = [19, 14]  # "R", "K"
        self.negativeAA = [15, 16]  # "D", "E"
        # "A", "I", "L", "M", "F", "W", "Y", "V"
        self.phobicAA = [8, 2, 1, 3, 0, 18, 9, 4]
        self.polarAA = [5, 7, 13, 12]  # "S", "T", "N", "Q"
        self.maxToxinLength = 3000  # default ???
        self.minToxinLength = 198  # I set this to the min len value of a single catalytic site
        self.maxAntidoteLength = 3000  # default ???
        self.minAntidoteLength = 198  # right now this minimum is arbitrary
        self.toxinSchema = ""  # the DNA sequence
        self.toxinTranslation = ""  # protein sequence
        # the code is translating the DNA and converting strings to an array of numbers
        self.toxinInteger = []
        self.antidoteSchema = ""
        self.antidoteTranslation = ""
        # the code is translating the DNA and converting strings to an array of numbers
        self.antidoteInteger = []
        # self.constraintLength = 50 #you can't allow infinitely large schemas
        if mutationRate == 0:
            self.mutationRate = random.randint(1, 30)
        else:
            self.mutationRate = mutationRate
        # tracking variables for catalytic states
        self.bestCatalyticDUBScore = 0
        self.bestCatalyticNucScore = 0
        self.bestNLSScore = 0
        self.bestTypeIVSSScore = 0
        # 0 is antidote and 1 is toxin, Just leave set to 0.5 to start
        self.NLSSiteLocation = 0.5
        # 0 is antidote and 1 is toxin, Just leave set to 0.5 to start
        self.TypeIVSiteLocation = 0.5
        self.combinedTAFitness = 0
        # setting up methods
        self.setSchemata(self.instantiationToggle)
        self.translateSchemata()
        self.translateNumbers()
        self.evaluateLocalBindingFitness()
        self.pairwiseAlignment(0)  # evaluate the DUB with pairwise
        self.pairwiseAlignment(1)  # evaluate the NUC with pairwise
        self.pairwiseAlignment(3)  # evaluate the NLS with pairwise
        self.pairwiseAlignment(4)  # evaluate the TYPEIV with pairwise
        self.sumCombinedTAFitness()

    # CLASS METHODS

    # Setters
    def setMutationRate(self, mutation_rate):
        self.mutationRate = mutation_rate

    def setMaxToxinLength(self, max_toxin_length):
        self.maxToxinLength = max_toxin_length

    def setMinToxinLength(self, min_toxin_length):
        self.minToxinLength = min_toxin_length

    def setMaxAntidoteLength(self, max_antidote_length):
        self.maxAntidoteLength = max_antidote_length

    def setMinAntidoteLength(self, min_antidote_length):
        self.minAntidoteLength = min_antidote_length

    def setToxinSchema(self, toxin_schema):
        self.toxinSchema = toxin_schema
        self.translateSchemata()
        self.translateNumbers()

    def setAntidoteSchema(self, antidote_schema):
        self.antidoteSchema = antidote_schema
        self.translateSchemata()
        self.translateNumbers()

    def re_evaluateAllFitness(self):
        self.evaluateLocalBindingFitness()
        self.pairwiseAlignment(0)  # evaluate the DUB with pairwise
        self.pairwiseAlignment(1)  # evaluate the NUC with pairwise
        self.pairwiseAlignment(2)  # evaluate the NLS with pairwise
        self.pairwiseAlignment(3)  # evaluate the TYPEIVSS with pairwise
        self.sumCombinedTAFitness()

    # A method to instantiate new TA schematas.
    # (This seeds with ORFs absent stop codons.)
    # 0 = instantiates perfect toxins, 1 = instantiates random toxins, 2 = instantiates a population from one mutated toxin.
    def setSchemata(self, toggle):
        toxinLength = random.randint(
            self.minToxinLength, self.maxToxinLength - 1)
        antidoteLength = random.randint(
            self.minAntidoteLength, self.maxAntidoteLength - 1)
        newToxinSchema = "ATG"
        newAntidoteSchema = "ATG"

        if toggle == 0:
            # USE THIS TO SEED WITH PERFECT TOXINS
            if random.randint(0, 1) == 0:  # make a nuclease
                newToxinSchema = "ATGGATCTG"
                newToxinSchema = newToxinSchema + \
                    random.choice(self.codonsNoSTOP) + "CTGCTG" + \
                    random.choice(self.codonsNoSTOP) + "CGT"
                for _ in range(0, 10):
                    newToxinSchema = newToxinSchema + \
                        random.choice(self.codonsNoSTOP)
                newToxinSchema = newToxinSchema + "CCGATTATTATTGAACTGAAA"
                for _ in range(0, 21):
                    newToxinSchema = newToxinSchema + \
                        random.choice(self.codonsNoSTOP)
                newToxinSchema = newToxinSchema + "GATCTGGTGCTG"
                for _ in range(0, 10):
                    newToxinSchema = newToxinSchema + \
                        random.choice(self.codonsNoSTOP)

                newToxinSchema = newToxinSchema + "CCGATTGGCCTGGAACTGAAA"
            else:  # make a dub
                newToxinSchema = "ATGCATTGGGTGACCCTGGTGATT"
                for _ in range(0, 9):
                    newToxinSchema = newToxinSchema + \
                        random.choice(self.codonsNoSTOP)
                newToxinSchema = newToxinSchema + "TATTAT" + \
                    random.choice(self.codonsNoSTOP) + "GATAGCCTG"
                for _ in range(0, 8):
                    newToxinSchema = newToxinSchema + \
                        random.choice(self.codonsNoSTOP)
                newToxinSchema = newToxinSchema + "ATT" + random.choice(self.codonsNoSTOP) + random.choice(
                    self.codonsNoSTOP) + random.choice(self.codonsNoSTOP) + "CTG"
                for _ in range(0, 5):
                    newToxinSchema = newToxinSchema + \
                        random.choice(self.codonsNoSTOP)
                newToxinSchema = newToxinSchema + "GAT"
                for _ in range(0, 9):
                    newToxinSchema = newToxinSchema + \
                        random.choice(self.codonsNoSTOP)
                newToxinSchema = newToxinSchema + "CAGCAG" + random.choice(self.codonsNoSTOP) + "GATGGC" + random.choice(
                    self.codonsNoSTOP) + random.choice(self.codonsNoSTOP) + random.choice(self.codonsNoSTOP) + "TGCGGC"
                for _ in range(0, 4):
                    newToxinSchema = newToxinSchema + \
                        random.choice(self.codonsNoSTOP)
                newToxinSchema = newToxinSchema + "GAAAAC"
            self.setToxinSchema(newToxinSchema)

        elif toggle == 1:
            # THIS SECTION INSTANTIATES RANDOM TA's WITH ONE SEEDWORD
            # This is adding random codons.
            # modulo 3 is because of codon triplets
            for _ in range(0, (toxinLength + 1)//3):
                newCodon = random.choice(self.codons)
                newToxinSchema = newToxinSchema + newCodon

            # this adds 1 seed word per toxin at a random index
            # why 360? I never see ORFs longer than ~120 Amino Acids
            index = random.randint(0, 360)
            # this was when I was mixing in seed words
            finalToxinSchema = newToxinSchema[0: index] + random.choice(
                self.seedWords) + newToxinSchema[index: len(newToxinSchema)]
            self.setToxinSchema(finalToxinSchema)

        elif toggle == 2:
            # THIS SECTION INSTANTIATES A POPULATION OF TOXINS FROM 1 FOUNDER CND BY MUTATION
            # Make a CND
            newToxinSchema = "ATGGATCTG"
            newToxinSchema = newToxinSchema + \
                random.choice(self.codonsNoSTOP) + "CTGCTG" + \
                random.choice(self.codonsNoSTOP) + "CGT"
            for i in range(0, 10):
                newToxinSchema = newToxinSchema + \
                    random.choice(self.codonsNoSTOP)
            newToxinSchema = newToxinSchema + "CCGATTATTATTGAACTGAAA"
            for i in range(0, 21):
                newToxinSchema = newToxinSchema + \
                    random.choice(self.codonsNoSTOP)
            newToxinSchema = newToxinSchema + "GATCTGGTGCTG"
            for i in range(0, 10):
                newToxinSchema = newToxinSchema + \
                    random.choice(self.codonsNoSTOP)
            newToxinSchema = newToxinSchema + "ATGCATTGGGTGACCCTGGTGATT"
            for i in range(0, 9):
                newToxinSchema = newToxinSchema + \
                    random.choice(self.codonsNoSTOP)
            newToxinSchema = newToxinSchema + "TATTAT" + \
                random.choice(self.codonsNoSTOP) + "GATAGCCTG"
            for i in range(0, 8):
                newToxinSchema = newToxinSchema + \
                    random.choice(self.codonsNoSTOP)
            newToxinSchema = newToxinSchema + "ATT" + random.choice(self.codonsNoSTOP) + random.choice(
                self.codonsNoSTOP) + random.choice(self.codonsNoSTOP) + "CTG"
            for i in range(0, 5):
                newToxinSchema = newToxinSchema + \
                    random.choice(self.codonsNoSTOP)
            newToxinSchema = newToxinSchema + "GAT"
            for i in range(0, 9):
                newToxinSchema = newToxinSchema + \
                    random.choice(self.codonsNoSTOP)
            newToxinSchema = newToxinSchema + "CAGCAG" + random.choice(self.codonsNoSTOP) + "GATGGC" + random.choice(
                self.codonsNoSTOP) + random.choice(self.codonsNoSTOP) + random.choice(self.codonsNoSTOP) + "TGCGGC"
            for i in range(0, 4):
                newToxinSchema = newToxinSchema + \
                    random.choice(self.codonsNoSTOP)
            newToxinSchema = newToxinSchema + "GAAAAC"
            self.setToxinSchema(newToxinSchema)

        for integers in range(0, (antidoteLength + 1)//3):
            newCodon = random.choice(self.codons)
            newAntidoteSchema = newAntidoteSchema + newCodon
        self.setAntidoteSchema(newAntidoteSchema)

    # A method to translate DNA code (4 bases to amino acids) using E.coli codons and 1 frame
    def translateSchemata(self):
        currentCodonIndex = 0
        currentTranslation = ""
        while currentCodonIndex < len(self.toxinSchema):
            codon = self.toxinSchema[currentCodonIndex: currentCodonIndex + 3]
            aminoAcid = self.switch_tRNA(codon)
            if aminoAcid == "*":
                currentTranslation = currentTranslation + aminoAcid
                break
            currentTranslation = currentTranslation + aminoAcid
            currentCodonIndex += 3
        self.toxinTranslation = currentTranslation

        currentCodonIndex = 0
        currentTranslation = ""
        while currentCodonIndex < len(self.antidoteSchema):
            codon = self.antidoteSchema[currentCodonIndex: currentCodonIndex + 3]
            aminoAcid = self.switch_tRNA(codon)
            if aminoAcid == "*":
                currentTranslation = currentTranslation + aminoAcid
                break
            currentTranslation = currentTranslation + aminoAcid
            currentCodonIndex += 3
        self.antidoteTranslation = currentTranslation

    # A method to translate AA code to numbers to speed up computational comparisons
    def translateNumbers(self):
        currentCodonIndex = 0
        currentTranslation = []
        while currentCodonIndex < len(self.toxinTranslation):
            aminoAcid = self.integer_Parser(
                self.toxinTranslation[currentCodonIndex])
            currentTranslation.append(aminoAcid)
            currentCodonIndex += 1
        self.toxinInteger = currentTranslation

        currentCodonIndex = 0
        currentTranslation = []
        while currentCodonIndex < len(self.antidoteTranslation):
            aminoAcid = self.integer_Parser(
                self.antidoteTranslation[currentCodonIndex])
            currentTranslation.append(aminoAcid)
            currentCodonIndex += 1
        self.antidoteInteger = currentTranslation

    def integer_Parser(self, aminoAcid):
        switcher = {
            "F": 0,
            "L": 1,
            "I": 2,
            "M": 3,
            "V": 4,
            "S": 5,
            "P": 6,
            "T": 7,
            "A": 8,
            "Y": 9,
            "*": 10,
            "H": 11,
            "Q": 12,
            "N": 13,
            "K": 14,
            "D": 15,
            "E": 16,
            "C": 17,
            "W": 18,
            "R": 19,
            "G": 20,
            "#": 21,
        }
        return switcher.get(aminoAcid)

    def switch_tRNA(self, codon):
        switcher = {
            "TTT": "F",
            "TTC": "F",
            "TTA": "L",
            "TTG": "L",
            "CTT": "L",
            "CTC": "L",
            "CTA": "L",
            "CTG": "L",
            "ATT": "I",
            "ATC": "I",
            "ATA": "I",
            "ATG": "M",
            "GTT": "V",
            "GTC": "V",
            "GTA": "V",
            "GTG": "V",

            "TCT": "S",
            "TCC": "S",
            "TCA": "S",
            "TCG": "S",
            "CCT": "P",
            "CCC": "P",
            "CCA": "P",
            "CCG": "P",
            "ACT": "T",
            "ACC": "T",
            "ACA": "T",
            "ACG": "T",
            "GCT": "A",
            "GCC": "A",
            "GCA": "A",
            "GCG": "A",

            "TAT": "Y",
            "TAC": "Y",
            "TAA": "*",
            "TAG": "*",
            "CAT": "H",
            "CAC": "H",
            "CAA": "Q",
            "CAG": "Q",
            "AAT": "N",
            "AAC": "N",
            "AAA": "K",
            "AAG": "K",
            "GAT": "D",
            "GAC": "D",
            "GAA": "E",
            "GAG": "E",

            "TGT": "C",
            "TGC": "C",
            "TGA": "*",
            "TGG": "W",
            "CGT": "R",
            "CGC": "R",
            "CGA": "R",
            "CGG": "R",
            "AGT": "S",
            "AGC": "S",
            "AGA": "R",
            "AGG": "R",
            "GGT": "G",
            "GGC": "G",
            "GGA": "G",
            "GGG": "G",
        }
        return switcher.get(codon, "*")

    def penaltyBindingFunction(self, currentMatches):
        # adding here a section that sets a maximum score for binding = 11. 11 is the sites reported in the crystal structure paper
        if currentMatches >= 11:
            return 1  # 1 if max binding is acheived
        else:
            return currentMatches/11

    def sumCombinedTAFitness(self):
        totalLength = len(self.toxinSchema) + len(self.antidoteSchema)
        if totalLength < 4500:  # magic coefficient loosely based on avg size of cif operons, but Parsimony pressure is on
            penalty = 0
        else:
            # magic coefficient loosely based on avg size of cif operons, but Parsimony pressure is on
            penalty = totalLength/4500
        self.combinedTAFitness = (self.bestCatalyticNucScore + self.bestCatalyticDUBScore +
                                  self.localBindingFitness + self.bestNLSScore + self.bestTypeIVSSScore - penalty)

    # The algorithm to use here is a simple sliding window.
    def evaluateLocalBindingFitness(self):
        # check for premature stop codon
        if len(self.toxinInteger) < 9 or len(self.antidoteInteger) < 9:
            self.localBindingFitness = -1000000
            return
        # Find the largest Schema
        if len(self.toxinInteger) >= len(self.antidoteInteger):
            largestTranslation = self.toxinInteger
            smallerTranslation = self.antidoteInteger
        else:
            largestTranslation = self.antidoteInteger
            smallerTranslation = self.toxinInteger

        # set point trackers
        maxFitness = -100  # arbitrary super low number
        numberSlides = 0
        movingTopPointer = 0
        movingBottomPointer = len(smallerTranslation) - numberSlides - 1
        while movingBottomPointer > -1:
            currentMatches = 0
            # compare numbers
            while movingBottomPointer < len(smallerTranslation):
                # check for salt bridges
                if largestTranslation[movingTopPointer] in self.positiveAA and smallerTranslation[movingBottomPointer] in self.negativeAA:
                    currentMatches += 1
                elif largestTranslation[movingTopPointer] in self.negativeAA and smallerTranslation[movingBottomPointer] in self.positiveAA:
                    currentMatches += 1

                # check for Van der Waals interaction
                # elif largestTranslation[movingTopPointer] in self.phobicAA and smallerTranslation[movingBottomPointer] in self.phobicAA:
                #    currentMatches += 2
                # check for hydrogen bonds (polar interaction)
                # elif largestTranslation[movingTopPointer] in self.polarAA and smallerTranslation[movingBottomPointer] in self.polarAA:
                #    currentMatches += 1

                # check for repulsion of charges
                elif largestTranslation[movingTopPointer] in self.positiveAA and smallerTranslation[movingBottomPointer] in self.positiveAA:
                    currentMatches -= 1
                elif largestTranslation[movingTopPointer] in self.negativeAA and smallerTranslation[movingBottomPointer] in self.negativeAA:
                    currentMatches -= 1

                # check for repulsion (hydrophilic with hydrophobic)
                # elif largestTranslation[movingTopPointer] in self.polarAA and smallerTranslation[movingBottomPointer] in self.phobicAA:
                #    currentMatches -= 1
                # elif largestTranslation[movingTopPointer] in self.phobicAA and smallerTranslation[movingBottomPointer] in self.polarAA:
                #    currentMatches -= 1

                movingTopPointer += 1
                movingBottomPointer += 1
                penaltyBindingFitness = self.penaltyBindingFunction(
                    currentMatches)
                if penaltyBindingFitness > maxFitness:
                    maxFitness = penaltyBindingFitness
                    self.localBindingFitness = penaltyBindingFitness
            numberSlides += 1
            # resets
            movingTopPointer = 0
            movingBottomPointer = len(smallerTranslation) - numberSlides - 1

        # reset pointers for other half of comparisons
        numberSlides = 1
        movingTopPointer = 0 + numberSlides
        movingBottomPointer = 0
        while movingTopPointer < len(largestTranslation):
            currentMatches = 0
            # compare numbers
            # problem here
            while movingBottomPointer < len(smallerTranslation) and movingTopPointer < len(largestTranslation):
                # check for salt bridges
                if largestTranslation[movingTopPointer] in self.positiveAA and smallerTranslation[movingBottomPointer] in self.negativeAA:
                    currentMatches += 1
                elif largestTranslation[movingTopPointer] in self.negativeAA and smallerTranslation[movingBottomPointer] in self.positiveAA:
                    currentMatches += 1

                # check for Van der Waals interaction
                # elif largestTranslation[movingTopPointer] in self.phobicAA and smallerTranslation[movingBottomPointer] in self.phobicAA:
                #    currentMatches += 2
                # check for hydrogen bonds (polar interaction)
                # elif largestTranslation[movingTopPointer] in self.polarAA and smallerTranslation[movingBottomPointer] in self.polarAA:
                #    currentMatches += 1

                # check for repulsion of charges
                elif largestTranslation[movingTopPointer] in self.positiveAA and smallerTranslation[movingBottomPointer] in self.positiveAA:
                    currentMatches -= 1
                elif largestTranslation[movingTopPointer] in self.negativeAA and smallerTranslation[movingBottomPointer] in self.negativeAA:
                    currentMatches -= 1

                # check for repulsion (hydrophilic with hydrophobic)
                # elif largestTranslation[movingTopPointer] in self.polarAA and smallerTranslation[movingBottomPointer] in self.phobicAA:
                #    currentMatches -= 1
                # elif largestTranslation[movingTopPointer] in self.phobicAA and smallerTranslation[movingBottomPointer] in self.polarAA:
                #    currentMatches -= 1

                movingTopPointer += 1
                movingBottomPointer += 1
                penaltyBindingFitness = self.penaltyBindingFunction(
                    currentMatches)
                if penaltyBindingFitness > maxFitness:
                    maxFitness = penaltyBindingFitness
                    self.localBindingFitness = penaltyBindingFitness
            numberSlides += 1
            # resets
            movingTopPointer = 0 + numberSlides
            movingBottomPointer = 0

    def toString(self):
        string = str(self.toxinTranslation) + "\n" + str(self.antidoteTranslation) + "\n" + "DUB fitness: " + str(self.bestCatalyticDUBScore) + "\n" + "NUC fitness: " + str(self.bestCatalyticNucScore) + "\n" + "NLS fitness: " + str(self.bestNLSScore) + "\n" + "TYPEIV fitness: " + \
            str(self.bestTypeIVSSScore) + "\n" + "Binding fitness: " + str(self.localBindingFitness) + "\n" + "combined fitness: " + \
            str(self.combinedTAFitness) + "\n" + "NLS Site: " + str(self.NLSSiteLocation) + \
            "\n" + "TYPEIV Site: " + str(self.TypeIVSiteLocation) + "\n"
        # print("universal toxin fitness: " + str(self.universalToxinFitness))
        # print("universal antidote fitness: " + str(self.universalAntidoteFitness))
        # print("universal fitness: " + str(self.universalFitness))
        print(string)
        return string

    # toggle is 0 for DUB, or 1 for NUC, 2 for NLS, 3 for TYPEIVSS
    def pairwiseAlignment(self, toggle):
        # Creating sample sequences
        # ref for this consensus of nuc2 site is Gillespie et al., Tangle Web
        DUB = Seq(
            "HWVTLVI---------YY-DSL--------I---L-----D---------QQ-DG---CG----EN")
        # ref for this consensus of nuc2 site is Gillespie et al., Tangle Web
        NUC = Seq(
            "DL-LL-R----------PIIIELK---------------------DLVL----------PIGLELK")
        NLS = "KRAR"  # from rossi et al., 1993 in agro Vird2
        # R-X(7)-R-X-R-X-R consensus from PNAS vergunst 2004
        TYPEIVSS = "R-------R-R-R"
        # checks the DUB score in toxin
        if toggle == 0:
            # divided by 23 to make perfect score 1 #this scoring sets opening gaps on the template to very minimal punishment
            self.bestCatalyticDUBScore = pairwise2.align.globalxs(
                DUB, self.toxinTranslation, -1, -.1, penalize_end_gaps=False, score_only=True)/23
        # checks the NUC score in toxin
        if toggle == 1:
            self.bestCatalyticNucScore = pairwise2.align.globalxs(
                NUC, self.toxinTranslation, -1, -.1, penalize_end_gaps=False, score_only=True)/23  # divided by 23 to make perfect score 1
        # checks the NLS score in both Toxin and Antidote
        if toggle == 2:
            self.bestNLSScore = pairwise2.align.globalxs(
                NLS, self.toxinTranslation, -1, -.1, penalize_end_gaps=False, score_only=True)/4  # divided by 4 to make perfect score 1
            self.NLSSiteLocation = 1  # 0 is antidote and 1 is toxin
            # now check the antidote
            # divided by 4 to make perfect score 1
            holderScore = pairwise2.align.globalxs(
                NLS, self.antidoteTranslation, -1, -.1, penalize_end_gaps=False, score_only=True) / 4
            if holderScore >= self.bestNLSScore:
                self.bestNLSScore = holderScore
                self.NLSSiteLocation = 0  # 0 is antidote and 1 is toxin
        # checks the TYPEIVSS score in both Toxin and Antidote
        if toggle == 3:
            self.bestTypeIVSSScore = pairwise2.align.globalxs(
                TYPEIVSS, self.toxinTranslation, -1, -.1, penalize_end_gaps=False, score_only=True)/4  # divided by 4 to make perfect score 1
            self.TypeIVSiteLocation = 1  # 0 is antidote and 1 is toxin
            # now check the antidote
            holderScore = pairwise2.align.globalxs(TYPEIVSS, self.antidoteTranslation, -1, -.1,
                                                   penalize_end_gaps=False, score_only=True)/4  # divided by 4 to make perfect score 1
            if holderScore >= self.bestTypeIVSSScore:
                self.bestTypeIVSSScore = holderScore
                self.TypeIVSiteLocation = 0  # 0 is antidote and 1 is toxin

Classes
    Main (this is just a main driver)
    TA   (this is the base of the model)
    CIEA (this runs the simulations)

Class Methods
    Main (Class)
        this just very simply instantiates a new CIEA class simulations and calls the runClassicTA() simulation method
    TA (Class)
        This class intantiates TA pairs as strings and converts those strings to integer codes.
        It contains all the evaluation methods.
    CIEA (Class)

        ///determine which parameters impact where NLS evolves.
        ///I should add some visualizer of where the domains are falling.
        ///just dubs or nucleases

        Priority
                ////REQUIRED FIXES 4////
        Finalize the logging functions.
        ///Fuck... I need to add functions to track average NLS fitness and average T4SS fitness



Version History:
1.0
    Was a numerical abstraction and was drastically re-implemented as actual DNA codons
1.1
    Introduced ATGC DNA codons and amino acid code translations instead of numerical code abstractions.
    Introduced graded charge interactions, hydrogen bonds, Van der Waals, polar/hydrophobic repulsion, ionic repulsion in the binding evaluation.
    Introduced a mutation rate parameter.
    Introduced an evolvable mutation rate.
    Introduced dub and nuc fitness calculations.
    Introduced minimal results logging functions.
    Introduced integer code conversions of strings in hope of speeding up algorithm.
    Tested parsimony pressure.
    Tested and Introduced parsimony pressure to DNA, but ORF pressure to translations.
    Tested feature to seed TAs with open reading frames.
        When I initialized with open reading frames it bogs the program down.
    Anecdotal Observations:
    Found that weight of parsimony pressure is having huge effects
    Found that to get the toxin sites, selection has to be stronger on the schema of the warhead than it does on binding.
        More specifically if selection is more or equally favoring binding, which is easier to evolve because the schema isn't defined per se, the warhead cannot remain intact.
    Found that when there is parsimony pressure the toxin and antidote reduces to nothing but what is necessary
    Found that if you shut off parsimony pressure - the binding region just keeps expanding and antidotes keep expanding.
        Without parsimony pressure is the only time I see length of toxin get longer than antidote.
    Found that I'm not maintaining diversity well. A single one is taking over the population before full toxins can evolve - one idea is to increase weight of toxin domain value in fitness - another is to up mutation rate. Also increasing mating pool and thereby decreasing killing selection.
    *this was fixed in later rounds with parameter optimization and recombination method tweaking.
1.2
    I will attempt to finalize and perfect the model in class (Research Methods in Evolutionary Computing, by Dr. Tauritz).
    Introduced JSON config file on 3.7.2022
    Introduced Random Seed defining for reproducibility on 3.7.2022
    I clocked the methods and the slowest methods were Mutation (taking ~3 seconds) and recombination (taking ~1 second).
        The way this works is recombination happens, then passes to mutation. Mutation can happen up to 30 times per individual (which slows) the process. Then mutation re-evaluates the fitness of all children which is slow.
        The conclusion is that fitness evals are taking the most time.
            The part of fitness that is taking the most time is binding evals, then DUB evals, then NUC evals.
                To speed up binding, test only charged residue analysis. This reduced the cost by about 1/3.
                    But it can never be better than N^2, because the inherent complexity is N2. This is probably as good as it gets.
                        I might be able to make it more efficient by including a max binding sites number concordant with the crystal structures (so the binding sites don't keep expanding).
    The default is now parsimony pressure and max binding score of 11
    Introduced constraint satisfaction with bonuses for subwords in catalytic domains. (this was later removed due to pairwise alignment)
    Introduced pairwise2 based alignments from biopython for all consensus domain predicting. This sped up evolution.
        The settings are : x- No parameters. Identical characters have score of 1, otherwise 0. this allows to not penalize for don't cares
                           s- The sequences have same open and extend gap penalties
    Introduced an instantiation Toggle where 0 starts simulation with perfect toxins and 1 starts simulation with de-novo random toxins
    Introduced double and triple point recombination.
    Introduced a termination condition based on average fitness not changing for a set number of generations.

1.3
    #make a results parser.
    test 0 parameter optimizations - DONE
    test 1 Instantiate de-novo
    test 2 Instantiate perfect toxins
    Test 3 Instantiate just Dub or Nuc
    Test 4 should encode actual real life CI toxins as instatiates, rather than evolving them against a fake template.
    Test 5 Test alternate sequences and sites
    Test 6 Figure out what factors change the sites
    Writeup
#!/usr/bin/env python3

import os
import sys
import random
import argparse
import copy
import glob


class Allele:
    '''
    Parent class of STRs, SNPs and indels
    subclasses are used to have custom allele types and mutation models
    the mutate method remains the same, however!
    '''
    
    def __init__(self, allele, probs):
       
        self.allele=allele
        self.numMuts=0
        self.origAllele=allele
        
        self.alphabet=[] # needs updating by other routines
        self.probs=[] # parallel array
        
        self.initAlleles(probs)      
 
    # both eq and hash assume that you are comparing/hashing within a locus.
    def __eq__(self, other):   
        return str(self.allele) == str(other.allele)

    def __hash__(self):
        return hash(self.allele)      
    
    def __str__(self):
        return self.allele

    def __copy__(self):
        '''
        copies, but copies deeply. (deepcopy's excessive bookkeeping is not needed. for our purposes, just use copy).
        Modified from:
        #https://stackoverflow.com/questions/1500718/how-to-override-the-copy-deepcopy-operations-for-a-python-object
        '''
        cls = self.__class__
        result = cls.__new__(cls)
        result.__dict__.update(self.__dict__)        
        result.alphabet = self.alphabet.copy()
        result.probs = self.probs.copy()
        return result
        
    def initAlleles(self, prob):
        '''
        This routine sets/updates the alphabet (and corresponding probabilities)
        Each kind of allele (snp, indel, str) has its own alphabet.
        '''
        pass
    
    def sequence(self, seqError):
        '''
        does not change state. returns an *estimate* of the true DNA sequence (str), given some (very small) sequencing error probability
        '''
        #alph=['A', 'C', 'G', 'T'] # (DNA alphabet lexically ordered; helps for evaluating the sequencing error routine)
        
        nerrors=0
        
        seqd=[]
        
        for c in self.allele:
            if random.random() < seqError:
                nerrors +=1
                # loop-less random character generator
                r = random.random()
                if r < 1/3.: # first letter in the alphabet (of the available to choose from)
                    seqd.append('C' if c == 'A' else 'A')
                elif r >= 2/3.: # last letter in the alphabet (of the available to choose from)
                    seqd.append('G' if c == 'T' else 'T')                    
                elif ord(c) <= ord('C'): # trickier case; C or A; the middle letter in the alphabet (subtracting out the true allele) is G
                    seqd.append('G')
                else: 
                    seqd.append('C')
            else:
                seqd.append(c)
        
        if nerrors:
            return "".join(seqd)
            
        return self.allele
        
        
    def mutate(self):
        '''
        all alleles can mutate; into what depends on the child class
        this is applied in PCR; ie, when copying
        note that mutate changes the allele state
        and the probability of mutation is taken to be constant
        '''
        rando = random.random()        
        newAllele=None
        cumulative = 0.
        i=0
        for nuc in self.alphabet:
            if nuc != self.allele: # note that the probabilities do not need to sum to 1 and you can omit the probability associated with the present allele. (makes updates easier)
                cumulative += self.probs[i]
                if cumulative > rando:
                    newAllele=nuc
                    self.numMuts +=1
                    break
            i += 1
                
        if newAllele is not None:
            self.allele=newAllele
            self.initAlleles(None) # None signals to not update the allele probabilties
            return 1
        return 0
        
    
class Snp(Allele):

    def initAlleles(self, probs): # probs is a float

        if probs is not None:
            self.alphabet = ['A', 'C', 'G', 'T'] # defaults for SNPs
            self.probs = [probs/3.] * 4 # probs is of type float; per round of PCR probability of transitioning to a new allele
        
        

class Indel(Allele):
    '''
    This is the most annoying class; the "alphabet" of indels is infinite.
    We take it as the set of characters seen in the population data
    (which you don't know until after you parse things)
    thus the flow is to make the different kinds of alleles
    then you update the indels so that they know the possible allelic states
    '''

    def updateAlphabet(self, alpha):
        self.alphabet = alpha       
        self.probs = [self.probs[0]] * len(self.alphabet)
        
    def initAlleles(self, probs): # probs is a float
    
        if probs is not None:
            self.probs = [probs] * max(1, len(self.alphabet)) # probs is of type float; per round of PCR probability of transitioning to a new allele           
        elif self.probs is None:
            print("indel error: Should not happen!\n", self, file=sys.stderr)
            exit(1)
           
class Str(Allele):

    def initAlleles(self, probs): # probs is a list of length 3
        assert probs is None or len(probs)==3
        
        # this determines the LUS sequence and its length
        offset = 4 # assumes tetramer. We may have to update 
        stop = len(self.allele) - offset + 1
        bestLen=0
        substrPos=0
                
        for shift in range(0, (offset-1)): # possible substring start locations of 0,1,2
            thisLen=1 # length of current match
            
            prevss= self.allele[shift:(offset+shift)] # first 4-mer
            i=0 # scope.
            
            for i in range((offset+shift), stop, offset):
                ss = self.allele[i:(i+offset)] # grab the next 4-mer 
                if ss == prevss: # same as last; it's bigger.
                    thisLen +=1

                elif bestLen < thisLen: # different from last and better; 
                    bestLen=thisLen
                    substrPos=i-offset
                    
                    thisLen=1
                else: # different from last and worse
                    thisLen=1
                    
                prevss=ss
            # trailing case; STR allele ends in the LUS
            
            if bestLen < thisLen:
               bestLen=thisLen
               substrPos=i-offset
        
        # motif associated w/ the LUS
        motif = self.allele[substrPos:(substrPos+offset)]
        
        # generate stutters:
                
        minusOne = self.allele[0:substrPos] + self.allele[ (substrPos+offset):]
        plusOne = self.allele[0:substrPos] + motif + motif + self.allele[ (substrPos+offset):]
        
        self.alphabet = [minusOne, self.allele, plusOne]

        if probs is not None:
            self.probs = probs # probs is of type float; per round of PCR probability of transitioning to a new allele

       
 
class MHap:
    '''
    Simple abstraction of a macrohaplotype
    as an ordered list of SNPs, indels, and STRs
    '''
    
    def __init__(self, locusString, locusName, snpMutProb, indelMutProb, minus1stutterProb, plus1stutterProb):
    
        self.locusString=locusString
        self.name=locusName
        self.initRecs(snpMutProb, indelMutProb, minus1stutterProb, plus1stutterProb)
        
    def __eq__(self, other):  
        if type(other)==str:
            return str(self.locusString) == other
            
        return self.locusString == other.locusString

    def __hash__(self):
        """Overrides the default implementation"""
        return hash(self.locusString)        
    
    def __str__(self):
        return self.locusString
        
    def __len__(self):
        L= len(self.snps)
        
        for i in self.indels:
            L += len(i)
        
        for i in self.strs:
            L += len(i)
        
        return L
        
    def __copy__(self):
        '''
        This is an attempt to speed things up; let's overwrite the copy.copy callback
        so that it is specific to this macrohaps (and not general)
        hopefully, in so doing, we can use copy.copy (instead of copy.deepcopy) and be both faster and accurate
        '''
        cls = self.__class__
        result = cls.__new__(cls)
        result.name = self.name
        result.locusString = self.locusString
        result.snps = [None] * len(self.snps)
        result.indels = [None] * len(self.indels)
        result.strs = [None] * len(self.strs)
                
        for i in range(len(self.snps)):
            result.snps[i] = copy.copy(self.snps[i])
            
        for i in range(len(self.indels)):
            result.indels[i] = copy.copy(self.indels[i])    
      
        for i in range(len(self.strs)):
            result.strs[i] = copy.copy(self.strs[i])           
            
        return result
        
    def initRecs(self,snpMutProb, indelMutProb, minus1stutterProb, plus1stutterProb):
        lstring = self.locusString
        # From Xuewen's code, a macrohaplotype is ;-separated, (snps;indels;strs), individual records within classes are csvs
        recs = lstring.split(";")
        if len(recs) != 3:
            print("Parsing problem with ", lstring, "too few record-types found!", file=sys.stderr)
            exit(1)
            
        #self.snps = recs[0].split(",")
        #self.indels = recs[1].split(",")
        #self.strs = recs[2].split(",")
        
        
        recs[0].split(",")
        self.snps= []
        for s in recs[0].split(","):
            snp = Snp(s, snpMutProb)
            self.snps.append(snp)
        
        self.indels=[]
        for i in recs[1].split(","):
            indel = Indel(i, indelMutProb)
            self.indels.append(indel)
            
        self.strs=[]
        for s in recs[2].split(","):
            STR = Str(s, [minus1stutterProb, 0., plus1stutterProb])
            self.strs.append(STR)
    
    def updateIndelI(self, indelIndex, alphabet):
        '''
        Indels don't know their alphabet; we learn that from data, which means we load the macrohap, then
        we update the indels so that each indel knows what states it's allowed to mutate into.
        This updates the alphabet for the ith indel; it'll update the mutation probabilities (uniform) too
        '''
        assert indelIndex < len(self.indels) 
        
        self.indels[indelIndex].updateAlphabet(alphabet)
    
    def sequence(self, sequencingError=0.0001):
        '''
        returns a string (in the locusString format)
        does NOT change state.
        this converts the true sequence into an estimated sequence 
        (a la MPS)
        a constant error probability across the DNA alphabet is assumed
        '''
        
        snps = [s.sequence(sequencingError) for s in self.snps]
        indels = [s.sequence(sequencingError) for s in self.indels]
        strs = [s.sequence(sequencingError) for s in self.strs]
        
        return ",".join(snps) + ";" + ",".join(indels) + ";" + ",".join(strs)
        
    
    def mutate(self):
        '''
        walk along the macrohap, and mutate each allele; if they mutate, you need to update the locus string
        returns the number of mutation events
        (note; changes state)
        ''' 
        
        nmuts=0
        for snp in self.snps:
            nmuts += snp.mutate()
            
            
        for indel in self.indels:
            nmuts += indel.mutate()

        for STR in self.strs:
            nmuts += STR.mutate()
        
        
        if nmuts >0:
            self.locusString = ",".join(str(s) for s in self.snps) + ";" + ",".join(str(s) for s in self.indels) + ";" + ",".join(str(s) for s in self.strs)
            
        return nmuts
        

def parseDB(f, snpError, indelError, minus1stutt, plus1stutt):
    
    if os.path.isfile(f):
        fh = open (f)
    elif f == '-':
        fh = sys.stdin
    else:
        print("Failed to open database: " , f , file=sys.stderr)
        exit(1)
    
    
    # population -> locus -> [[allele, count]]
    d = {}
    
    indelStates={}
    j=0
    for line in fh:
        sp = line.rstrip().split("\t")
        if len(sp) < 4:
            continue
        loc = sp[0] # locus id
        pop = sp[1] # the population
        
        if pop not in d:
             d[pop] = {}
            
        outer = d[pop]
        
        
        
        if loc not in outer:
            outer[loc] = []
            indelStates[loc] = []
            
        mhap = MHap(sp[2], loc, snpError, indelError, minus1stutt, plus1stutt)

        outer[loc].append( [mhap, int(sp[3]) ] )
        
        j=0
        locIndels = indelStates[loc]
        
        for indel in mhap.indels:
            if len(locIndels)==j:
                locIndels.append( set() )
                
            locIndels[j].add( indel.allele) 
            j += 1
        
    for loc, states in indelStates.items():
        j=0
        for j in range(len(states)): # for all j indels
            #print(loc, state, j)
            states[j]=list(states[j])
            j +=1
            
    # I need to update the feasible alleles for each indel ineach macrohap... lots'o'for!
    for pop, outer in d.items():
        for loc, mhaps in outer.items():
            for mhap, count in mhaps:
                j=0
                for indel in mhap.indels:
                    mhap.updateIndelI(j, indelStates[loc][j])
                    j += 1
    
    if fh != sys.stdin:
        fh.close()

    return d
   

def getRandomHap(macrohaps):
    '''
    takes in the 2D list (part of the macrohaplotype database)
    and samples a haplotype at random
    '''
    cumulative =0
    for hap, count in macrohaps:
        cumulative += count
    
    rando = random.randint(0, cumulative)
    cumulative =0
    for hap, count in macrohaps:
        cumulative += count
        if cumulative >= rando:
            return hap
    # cannot happen.
    return None
    
def applyPcr(alleleCounts, probCopied):
    '''
    Performs 1 round o PCR on a single locus
    takes in a dictionary (macrohap -> count)
    and returns a dictionary (macrohap -> count, after 1 round of PCR)
    PCR increases the number of alleles, and with mutation, the number of distinct haplotypes as well
    '''
    
    out={}
    nmuts =0
    for allele, count in alleleCounts.items():
        

        # this is a bit inefficient, but so be it. it's python anyways; you want efficiency, use another language!
        # note: copy has been overloaded to copy deeply (and efficiently).
        # this lowered the runtime (somewhat)
        clone = copy.copy(allele)    
        #clone = copy.deepcopy(allele)  
        #clone = pickle.loads(pickle.dumps(allele, -1))
        for i in range(count):
            
            if allele in out:
                out[allele] += 1    
            else:
                out[allele]=1
        
            if random.random() >= probCopied: # not copied in PCR
                continue # count remains the same

           
            # copied,
            if clone.mutate()>0: # and at least 1 mutation

                if clone in out:
                    out[clone]+=1
                else: 
                    out[clone]=1
                
                #clone = copy.deepcopy(allele)
                clone = copy.copy(allele) 
                #clone = pickle.loads(pickle.dumps(allele, -1))
            else: # copied w/o pcr error
                out[allele] += 1    # key must exist.


    return out


def sequencingSimulator(state, sequencingErrorProb, samplingProb):
    loci = list(state.keys())
    finalState = {}
      
    for loc in loci:
        out = {}
        finalState[loc] = out
        for allele, count in state[loc].items():
            for i in range(count):
                if random.random() < samplingProb: # was this particular PCR copy sequenced?
                    seq = allele.sequence(sequencingErrorProb) # and let's sprinkle in some sequence error. note that seq is a string.
                    if seq in out:
                        out[seq]+=1
                    else:
                        out[seq]=1
                        
    return finalState
    
def pcrSimulator(initialState, probCopied, nrounds=30):
    
    loci = list(initialState.keys())
    finalState = {}

    for r in range(nrounds):   
        #print("Round " , r, file=sys.stderr)
        
        for loc in loci:
            if loc not in finalState:
                finalState[loc] = applyPcr(initialState[loc], probCopied)
            else:
                finalState[loc] = applyPcr(finalState[loc], probCopied)
    
        
    
    
    return finalState


def getTrueGenotypes(db, genotypes=None, locusFilter=None):
    '''
    This generates a genotype assuming Hardy Weinberg (ie, including, no mutation, male/female alleles are independent)
    in practice, you randomly sample two alleles with replacement.
    This returns a dictionary of lists of lists;
    the key is the locus -> outer list
    there is one entry (for each locus) for each individual (outerlist)
    and each inner list is a pair of alleles (together, the genotype)
    
    This function can be run multiple times; each time it adds the new genotypes to the list
    '''
    if genotypes is None:
        genotypes = {}

    if locusFilter is not None:
        locusFilter = locusFilter.upper()
    
    nloc=0
    for locus, haps in db.items():
        h1 = getRandomHap(haps)
        h2 = getRandomHap(haps)
        
        if locusFilter is not None and locus.upper() != locusFilter:
            continue
            
        nloc +=1
        if locus not in genotypes:
            genotypes[locus] = [[h1,h2]]
        else:
            genotypes[locus].append([h1,h2])
        
    if nloc==0:
        if locusFilter is None:
            print("Some problem reading the allele database...", file=sys.stderr)
            exit(1)
            
        for locus, haps in db.items():
            print(locus, file=sys.stderr)
            
        print("No locus " , locusFilter, "found." , "The above loci are the valid loci to choose from", file=sys.stderr, sep="\n")
        exit(1)
        
    return genotypes


def printSim(sim, initialLoadings, argv, out):
    
    print("#Markername\tCounts\tHapVarLen\thapVar(s)\tLoading", file=out)
    
    for locus, inner in sim.items():
        for allele, count in sorted(inner.items(), key=lambda item: -item[1]):
            loading=0
            if allele in initialLoadings[locus]:
                loading = initialLoadings[locus][allele]
                del initialLoadings[locus][allele]
            # , is the SNP/indel/str delimiter  (within)
            # ; is the delimiter between ( always 2)
            print(locus, count, len(allele)-allele.count(",")-2, allele, loading, sep="\t", file=out)
        for allele, loading in initialLoadings[locus].items():
            s = str(allele)
            print(locus, 0, len(s)-s.count(",")-2, s, loading, sep="\t", file=out)

def loadNCells(basedir, args, cellLoadings):
    '''
    Reads from the cache;
    specifies the population (eg, CEU)
    and a vector of cell loadings (ints)
    eg, [2]
    would be a single source sample, 2 replicates; ie, 2 cells worth of DNA from the same person
    while 
    [2,5,6]
    would be a 3 person mixture with 2, 5 and 6 cells worth of DNA for each of the three contributors
    
    '''

    out = {}
    initialLoadings = {}
    
    pop = args.Pop.upper()

    if not  os.path.isdir( os.path.join( basedir, pop)):
        print("Failed to find directory : ", basedir , pop , sep="\n", file=sys.stderr)
        exit(1)

    # get all of the files...
    files = glob.glob( os.path.join(basedir, pop, "cell.*tsv"))
    # filter out the empty files (the cache is writing. technically this can fail if a file is half-written; let's ignore that for now
    files = [ f for f in files if os.stat(f).st_size>0]

    #  'cell.20.42.CEU.tsv'
    # is the format; means this is replicate 42 of cell 20.
    # a cell corresponds to some PCR outcome based on a single diploid cell's alleles
    # ...
    # we want to sample (without replacement) some number of replicates from some cell
    
    s = {}
    for f in files:
        cell = f.split(".")[-4] # see 5 lines up; file format;
        if cell not in s:
            s[cell]=[f]
        else:
            s[cell].append(f)
    
    # heavy handed, but ensure we have enough replicates to sample from
    x = max(cellLoadings)
    
    keys = list(s.keys())
    # we have a cell; make sure we have enough replicates to sample.
    for k in keys:
        if len(s[k]) < x:
            del s[k]
    
   
    if len(s.keys()) < len(cellLoadings):
        print("You need more simulations for:", cellLoadings, pop, sep="\n", file=sys.stderr)
        exit(1)
    
    # apply the second PRNG seed. this gives us random PCR outcomes (based onthe individual specified in the first seed value)
    if len(args.R)>1:
        random.seed(args.R[1])
    # these are the random individuals that we've sampled
    inds = random.sample( list(s.keys()), len(cellLoadings)) 
    
    # now we need to select which cells contributed from these individuals
    files = []
    
    i=0
    for ind in inds:
        #print(ind)
        files.extend( random.sample( s[ind], cellLoadings[i]))
        i+=1
   
    locFilt = None
    if args.L != None:
        locFilt = args.L.upper()
    
    for f in files:
        with open(f) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                s = line.rstrip().split("\t")
                locus = s[0]
                # this is where we constrain to a specified locus.
                if locFilt is not None and locus.upper() != locFilt:
                    continue
                obs = int(s[1])
                allele = s[3]
                
                loading = int(s[4])
                
                if locus not in out:
                    out[locus] = {}
                
                inner = out[locus]
                
                if allele not in inner: # the string version and the mhap object are equivalent in the dictionary; we want the keys to be mhap objects.
                    mhap = MHap(allele, locus, 0., 0., 0., 0.) # pcr mutation properties are set to 0; (that's because PCR has aleady happened, so why bother w/ what the values should be)
                    inner[mhap] = obs
                else: # but strings suffice for incrementing the value
                    inner[allele] += obs
                
                if loading>0: # this is the true number of molecules, prior to PCR. 
                
                    if locus not in initialLoadings:
                        initialLoadings[locus]={}
                    inner=initialLoadings[locus]
                    
                    if allele not in inner:
                        inner[allele]=loading
                    else:
                        inner[allele] += loading
                    
                
   
    return out, initialLoadings
    


def simulator(args):

    pg_per_diploid_genome = 6.6

    defaultdb = os.path.abspath(os.path.dirname(__file__))
   
    defaultout = defaultdb
    # these are the sites that are internally consistent between TRcaller (MiSeq; Forenseq) and the 1KG project (shotgun).
    defaultdb = os.path.join(defaultdb, "../../data/stitchedMhAlleles/stitched.MH.allele_count.1KG2504.pop.trcaller_consistent.tsv")
    # used when building caches...
    defaultout = os.path.join(defaultout, "../../data/simulations")
    
    parser = argparse.ArgumentParser(description="Let's make some simulations!")
    
    # PCR parameters
    parser.add_argument("-R", "--random_number_generator_seed", dest='R', nargs="+", help="Used to seed the PRNG; 1st arg: seeds the genotype; 2nd arg; seeds mutation process; both are optional", default=[], type=int)
 
    parser.add_argument("-r", "--rounds_of_pcr", dest='r', help="The number of rounds of PCR (default: 30)", default=30, type=int)
    parser.add_argument("-e", "--pcr_effiency", dest='e', help="PCR efficiency (ie, probability of copying), per round (default: 0.672242)", default=0.672242, type=float)
    parser.add_argument("-s", "--snp_pcr_error", dest='s', help="SNP/indel error rate, per PCR round  (default: 0.000018)", default=0.000018, type=float)
    parser.add_argument("-m", "--minus1_stutter_error", dest='m', help="minus 1 stutter rate; per PCR round (default: 0.0035)", default=0.003505857, type=float)
    # taken as: 1-((1-0.1)**(1/30))
    # where 30 is the number of PCR rounds (which is what we use), and a 10% stutter ratio is expected (this is a rough estimate that doesnt consider back mutation)   
    parser.add_argument("-p", "--plus1_stutter_error", dest='p', help="plus 1 stutter rate; per PCR round (default: 0.000334)", default=0.0003349551, type=float)
    # taken as: 1-((1-0.01)**(1/30))
 
    # sequencing information
    parser.add_argument("-Q", "--sequencing_error", dest='Q', help="Per nucleotide sequencing error rate (defaults to q30; 0.001)", default=0.001, type=float) # q40
    parser.add_argument("-E", "--sequencing_effiency", dest='E', help="What fraction of the generated sequences get sequenced? (defaults to 0.01)", default=0.01, type=float) 
    
    # locus info
    parser.add_argument("-L", "--locus", dest='L', help="Restrict simulation to locus L (defaults to all loci)", default=None, type=str) 
    
    
    # population information
    parser.add_argument("-D", "--population_database", dest='d', help="file of population counts...", default=defaultdb, type=str)   
    
    # sample information
    parser.add_argument("-P", "--population_to_simulate", dest='Pop', help="The population to simulate; use 1KG/Coriell 3-letter identifiers. Default: CEU", default='CEU', type=str)   
    parser.add_argument("-N", "--number_of_contributors", dest='N', help="The number of individuals simulated", default=1, type=int)
    parser.add_argument("-S", "--starting_DNA_pg", dest='S', help="The number of picograms of DNA", default=pg_per_diploid_genome, type=float)
    parser.add_argument("-C", "--starting_DNA_cells", dest='C', help="The number of cells of DNA", default=0, type=int)
    parser.add_argument("-F", "--mixture_fractions", dest='F', nargs="+", help="The mixture fractions; one for each N", default=[], type=float)
 
    # output information
    parser.add_argument("-B", "--build_cache", dest='B', help="Creates a cache...", action='store_true', default=False)   
    parser.add_argument("-U", "--use_cache", dest='U', help="Uses cached data...", action='store_true', default=False) 
 
    (results, args) = parser.parse_known_args(args[1:])[0:2]
    
    if len(args):
        print("Extra arguments detected...", args, sep="\n", file=sys.stderr)
        return 1

    outstream = sys.stdout
    if results.B:
        # caches simulate a single cell's worth of DNA 
        # with 100% sequencing efficiency
        results.C=1
        results.N=1
        results.E=1.0
        results.F =[1.0]
        results.Q=0.0 # turn off sequencing error. it's only used when you sequence from the cache
        
        # need determinism
        if len(results.R)==0:
            results.R.append( random.randint(0, 1e6))
        if len(results.R)==1:
            results.R.append( random.randint(0, 1e6))
        
        defaultout = os.path.join(defaultout, results.Pop)
        
        
        if not os.path.isdir( defaultout):
            print(defaultout , "is not a valid output directory. Please make it!", file=sys.stderr, sep="\n")
            return 1
        
        outfile = os.path.join(defaultout, "cell." + ".".join([str(seed) for seed in results.R]) + "." + results.Pop + ".tsv")
        if os.path.isfile(outfile):
            print(outfile , "already exists. I suggest using different seeds!", file=sys.stderr, sep="\n")
            return 1
            
        outstream = open(outfile, "w")
        print("#" , results , file=outstream)
        return 0
        
    fracs = results.F
    
    # assume equal parts mixtures, if not specified
    if len(fracs)==0:
        fracs = [1.0/results.N]*results.N
    elif len(fracs) > results.N:
        print("Too many mixture fractions requested!", file=sys.stderr)
        return 1
    elif len(fracs) < results.N-1:
        print("Too few mixture fractions requested!", file=sys.stderr)
        return 2
        
    cumu = 0
    for f in fracs:
        cumu +=f
        if f < 0:
            print("negative mixture fraction?!", file=sys.stderr)
            return 3
            
    # can specify N-1 mixture fractions...
    if len(fracs) == results.N-1:
        fracs.append(max(0, 1-cumu))
    
    # taken as a real number; how many haploid genomes (in total) are we sampling
    if results.C>0:
        totalLoading = results.C # NOT times two as it's applied once per haploid genome (below) (for gt in gtpair)
    else:
        totalLoading = (results.S / pg_per_diploid_genome)
        
    if min([totalLoading*f for f in fracs]) < 0.9: # allow a modest amount of rounding error. purest would say it should be 1.
        print("Your mixture hypothesis is incompatible with the number of cells sampled... that's not good!", file=sys.stderr)
        return 1

    if len(results.R)>0:
        random.seed(results.R[0])
    
    if results.U: # routine for loading data from the cache...
        #print("fracs", fracs)
        #def loadNCells(basedir, args, cellLoadings):
        # note that the only substantive difference is that macrohap alleles cannot be PCR'd (their mutation rates are set to 0)
        sim, initialLoadings = loadNCells( defaultout, results, [round(f*totalLoading) for f in fracs]) 
               
    else:
        #parseDB(f, snpError, indelError, minus1stutt, plus1stutt, normalize=True):
        db = parseDB(results.d, results.s, results.s, results.m, results.p)
        
        pop = results.Pop.upper()
        
        if pop not in db:
            print("Failed to find pop: ", pop, "The following populations are available", sep="\n", file=sys.stderr)
            for pop in db.keys():
                print(pop, file=sys.stderr)
            return 1
            
        # locus -> [ [h1,h2],[h3,h4] ] (for a 2-person mixture)
        trueGts = {}
        for i in range(results.N):
            trueGts = getTrueGenotypes(db[pop], trueGts, results.L)
            
        # locus -> dictionary of alleles and their count
        initialLoadings={}
        for locus, gts in trueGts.items():
            initialLoadings[locus] = {}
            inner = initialLoadings[locus]
            i=0
            for gtpair in gts:
            
                for gt in gtpair:
                    # the user is given a lot of latitude. they can ask for abstractions that just make no sense. (fractional cells worth of DNA).
                    if gt not in inner:
                        inner[gt] = totalLoading * fracs[i]

                    else:
                        inner[gt] += totalLoading * fracs[i]
                        
                
            # we're counting absolute numbers of molecules. this must be discrete.
            for allele in inner:
                inner[allele] = max(1, round(inner[allele]))
                #print(str(allele), inner[allele])
        
        if len(results.R)>1:
            random.seed(results.R[1])
            
        sim = pcrSimulator(initialLoadings, results.e, results.r)
        
    #def sequencingSimulator(state, sequencingErrorProb, samplingProb):
    sim = sequencingSimulator(sim, results.Q, results.E)
    
    #randomly samples alleles (sequencing efficiency) and adds sequencing error
    printSim(sim, initialLoadings, args, out=outstream)
    
    if outstream != sys.stdout:
        outstream.close()
        
    
    return 0
    


if __name__ == "__main__":
  exit( simulator(sys.argv) )

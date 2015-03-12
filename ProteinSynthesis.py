import random
import socket
import sys
from Bio import Entrez
from Bio import SeqIO
from BioSQL import BioSeqDatabase
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord

#Structure of each cells functions, time and quantity
class Cell:

    def __init__(self, id, dnaSeq, maxPSNum, maxMutNum):
        self.id = id
        self.startDNASeq = dnaSeq
        self.currDNASeq = list(dnaSeq)
        self.mrnaSeq = ""
        # self.ribosome = Ribosome()
        # self.mRNA = mRNA()
        # self.tRNA = tRNA()
        self.atpAmount = 0
        self.proteins = []
        self.origReps = []
        self.errorRate = 0
        self.timeToLive = 110
        self.DNAState = 0
        self.timeToRep = 150
        self.replicate = False
        self.dnaDiff = 0
        self.ps_amtTime = 45
        self.ps_srtNum = random.randrange(0, maxPSNum+1)
        self.mt_Num = random.randrange(0, maxMutNum+1)
        self.numPS = 0
        self.numMuts = 0

    def Respiration(self):
        None

    def Glycolysis(self):
        None

    def Transcription(self):
        '''
        Coding Strand (5' -> 3'):
        ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG

        Template Strand (3' -> 5'):
        TACCGGTAACATTACCCGGCGACTTTCCCACGGGCTATC

        Transcription works on the Template Strand in the 5' to 3' direction
        by performing a reverse complement to get the mRNA (5' -> 3'):

        AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG

        In BioPython however we just work with the coding strand and get the mrna by
        switching all T -> U.
        '''

        coding_strand = Seq(''.join(self.currDNASeq), IUPAC.unambiguous_dna)
        template_strand = coding_strand.reverse_complement()

        # Bio Python method
        self.mrnaSeq = str(coding_strand.transcribe())

        # TRue Biological Process
        mRNA_strandB = template_strand.reverse_complement().transcribe()


    def Translation(self):
        # Find Proteins attributed to a mrnaSeq
        messenger_rna = Seq(str(self.mrnaSeq), IUPAC.unambiguous_rna)

        # Stop translation at first in frame stop codon as described in dna table
        messenger_rna.translate()

        # We can also specify our own stop codon
        # Use for our Simulation by splicing dna seq at various points?
        # messenger_rna.translate(stop_symbol="@")


    def ProteinSynthesis(self):

        self.rnaSeq = []
        self.numPS += 1

        # Transcribe dna strand and find complement rna sequence via RNA polymerase
        self.Transcription()

        self.Translation()

    # Define Death of a cell by (changing graph color and removing resources and setting weight to 0)
    def Apoptosis(self):
        None

    # Mutate dna sequence by lering nucleotides randomly
    def DNA_Damage(self):

        mutIndex = random.randrange(len(self.currDNASeq))

        if self.currDNASeq[mutIndex] == 'A':
            self.currDNASeq[mutIndex] = 'U'
        elif self.currDNASeq[mutIndex] == 'U':
            self.currDNASeq[mutIndex] = 'A'
        elif self.currDNASeq[mutIndex] == 'C':
            self.currDNASeq[mutIndex] = 'G'
        elif self.currDNASeq[mutIndex] == 'G':
            self.currDNASeq[mutIndex] = 'C'
        else:
            self.currDNASeq[mutIndex] = '_'

        self.numMuts += 1

    def DNA_Repair(self):
        None

    def __iter__(self):
        return self.__iter__()

def ProteinSynthesis_Sim(numSeconds, dnaSeq, maxPSNum, numCells, maxMutNum):

    ps_cellLists = []

    # Create a hash table by using a list to hold lists of cells indexed by their respective ps_startNum
    # Pro: Efficiently begin the PS process for all cells with the same startNum
    # Can use knight's tour to begin process for a group of cells with same startNum
    for i in range(maxPSNum+1):

        ps_cellLists.insert(i, [])

    for i in range(numCells):

        aCell = Cell(i, dnaSeq, maxPSNum, maxMutNum)
        print(aCell.ps_srtNum)
        ps_cellLists[aCell.ps_srtNum].append(aCell)

    for currentSecond in range(numSeconds):

        #  Check if a cell is ready to perform protein synthesis
        psTime = GenPSTime(maxPSNum)

        mutateTime = GenPSTime(maxPSNum)

        for i in range(psTime):
            for cell in ps_cellLists[i]:

                cell.ProteinSynthesis()

                #  Check if a cell had a mutation occur
                if cell.mt_Num == mutateTime:
                    cell.DNA_Damage()

    for i in range(maxPSNum):
        for cell in ps_cellLists[i]:
            print("Cell #" + str(cell.id) + "(Num PS: " + str(cell.numPS) + ", Num Damage: " + str(cell.numMuts) + "): " + str(cell.currDNASeq))

def GenPSTime(maxPSNum):

    return random.randrange(0, maxPSNum+1)



# Create A Menu to display
showMenu = True
menuSelection = 0

# myBioServer = BioServer()

# Create cells and run program to simulate
# program has pieces which move in various ways through the knights tour problem
# but are dictated by AI logic as to which path to take

ProteinSynthesis_Sim(100, "AGCT", 200, 5, 200)

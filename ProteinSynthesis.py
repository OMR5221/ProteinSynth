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

# Each Vertex uses a dictionary to hold adjacent vertices and weight
class Vertex:
    def __init__(self, key):
        self.id = key
        self.vertsConnectedTo = {}

    def addNeighbor(self, nbrVert, connectWeight):

        self.vertsConnectedTo[nbrVert] = connectWeight

    def __str__(self):
        return str(self.id) + " connected to: " + str([x.id for x in self.vertsConnectedTo])

    # Returns the id for all neighbors for this Vertex
    def getNeighbors(self):
        return self.vertsConnectedTo.keys()

    def getID(self):
        return self.id

    def getWeight(self, nbrVert):
        return self.vertsConnectedTo[nbrVert]

# ADT to hold Graph framework: Adjacency List (Used for coordinating knights tour)
# Holds a dictionary of all Vertices
class AL_Graph:

    def __init__(self):
        # Key: Vertex ID, Value: Vertex Object
        self.vertexList = {}
        self.numVertices = 0

    def __contains__(self, n):
        return n in self.vertexList

    # Make iterating through VertexList values possible
    def __iter__(self):
        return iter(self.vertexList.values())

    # Creates and adds a Vertex Node to Graph's Vertex List
    def addVertex(self, key):

        newVert = Vertex(key)

        # Add new Verte to this Graph
        self.vertexList[key] = newVert

        self.numVertices += 1

        return newVert

    # Find a particular vertex in the graph using key
    def getVertex(self, vertKey):

        if vertKey in self.vertexList:
            return self.vertexList[vertKey]

        else:
            return None

    # Define a weighted connection between two Vertices
    def addEdge(self, vert1Key, vert2Key, weight):

        if vert1Key not in self.vertexList:
            # Add the Verte to the Graph
            vert1 = self.addVertex(vert1Key)

        if vert2Key not in self.vertexList:
            vert2 = self.addVertex(vert2Key)

        # Update vert1.vertsConnectedTo dictionary
        self.vertexList[vert1Key].addNeighbor(self.vertexList[vert2Key], weight)

    # Returns a list of all Vertices by their key
    def getVertices(self):
        return self.vertexList.keys()

#Structure of each cells functions, time and quantity
class Cell:

    def __init__(self, id, dnaSeq, maxPSNum, maxMutNum):
        self.id = id
        self.startDNASeq = dnaSeq
        self.currDNASeq = dnaSeq
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
        self.ps_srtNum = random.randrange(1, maxPSNum+1)
        self.mt_Num = random.randrange(1, maxMutNum+1)
        self.numPS = 0

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

        coding_strand = Seq(self.currDNASeq, IUPAC.unambiguous_dna)
        template_strand = coding_strand.reverse_complement()

        # Bio Python method
        self.mrnaSeq = str(coding_strand.transcribe())

        # TRue Biological Process
        mRNA_strandB = template_strand.reverse_complement().transcribe()


    def Translation(self):
        # Find Proteins attributed to a mrnaSeq
        messenger_rna = Seq(str(self.mrnaSeq), IUPAC.unambiguous_rna)

        # Stop translation at first in frame stop codon as described in dna table
        messenger_rna.translate(to_stop=True)

        # We can also specify our own stop codon
        # Use for our Simulation by splicing dna seq at various points?
        messenger_rna.translate(stop_symbol="@")


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
    def DNA_Damage(self, mutateNum):

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

    def DNA_Repair(self):
        None

    def __iter__(self):
        return self.__iter__()

def ProteinSynthesis_Sim(numSeconds, dnaSeq, maxPSNum, numCells, maxMutNum):

    ps_cellLists = []

    # Create a hash table by using a list to hold lists of cells indexed by their respective ps_startNum
    # Pro: Efficiently begin the PS process for all cells with the same startNum
    # Can use knight's tour to begin process for a group of cells with same startNum
    for i in range(maxPSNum):

        ps_cellLists.insert(i, [])

    for i in range(numCells):

        aCell = Cell(i, dnaSeq, maxPSNum, maxMutNum)

        ps_cellLists[aCell.ps_srtNum].append(aCell)

    for currentSecond in range(numSeconds):

        #  Check if a cell is ready to perform protein synthesis
        psTime = GenPSTime(maxPSNum)

        mutateTime = GenPSTime(maxPSNum)

        for i in range(psTime):
            for cell in ps_cellLists[i]:
                if cell.ps_srtNum == psTime:
                    cell.ProteinSynthesis()
                #  Check if a cell had a mutation occur
                elif cell.mt_Num == mutateTime:
                    cell.DNA_Damage()

    for i in range(maxPSNum):
        for cell in ps_cellLists[i]:
            print("Cell " + str(cell.id) + "(" + str(cell.numPS) + "): " + str(cell.mrnaSeq))

def GenPSTime(maxPSNum):

    return random.randrange(1, maxPSNum+1)



# Create A Menu to display
showMenu = True
menuSelection = 0

# myBioServer = BioServer()

while showMenu:
    print("\n\n")
    print("Welcome to Protein Synthesis Simulation!\n")
    print("Please select an option below.\n\n")
    print("1. Connect to a BioSQL Server\n")
    print("2. Create a Sub Database\n")
    print("3. Connect to a Sub Database\n")
    print("4. Fetch and Load a Sequence\n")
    print("5. Run Protein Synthesis\n")
    print("\n")

    menuSelection=int(raw_input())

    if menuSelection == 1:
        print("Connecting to the specified Server: \n")
        ''' if myBioServer.connected == False:
            print("Initiating BioSQL Server Connection Process:\n")
            myBioServer.Connect()
        else:
            print("You are already connected to a BioSQL Server: \n")
            print("\t%s as user %s" % (myBioServer.serverName, myBioServer.user)) '''

    elif menuSelection == 2:
        print("Creating a Sub Database: \n")
        # myBioServer.Create_SubDB()

    elif menuSelection == 3:
        print("Connecting to Sub Database: \n")
        ''' myBioServer.Connect_Server()
        myBioServer.Connect_DB()'''

    elif menuSelection == 4:
        print("Fetching and Loading a Sequence: \n")
        '''myBioServer.Connect_Server()
        myBioServer.Connect_DB()
        myCell.currDNASeq = myBioServer.FetchSeq()'''

    elif menuSelection == 5:
        print("Starting Protein Synthesis Simulation: \n")
        ProteinSynthesis_Sim(100, "AGCT", 200, 15, 200)

    else:
        print("Error")


# Create cells and run program to simulate
# program has pieces which move in various ways through the knights tour problem
# but are dictated by AI logic as to which path to take
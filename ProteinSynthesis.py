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

        # Add new Vertex to this Graph
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
            # Add the Vertex to the Graph
            vert1 = self.addVertex(vert1Key)

        if vert2Key not in self.vertexList:
            vert2 = self.addVertex(vert2Key)

        # Update vert1.vertsConnectedTo dictionary
        self.vertexList[vert1Key].addNeighbor(self.vertexList[vert2Key], weight)

    # Returns a list of all Vertices by their key
    def getVertices(self):
        return self.vertexList.keys()

# Look at all potential dna sequences off of this cell's current one
# Figure out the smallest number of transformations needed to turn the starting word into the ending word
def seqLadder(dnaFile):

    dnaDict = {}

    dnaLadderGraph = AL_Graph()

    # Read in a sequence to generate all possible iterations of the sequence
    dnaFileIN = open(dnaFile, 'r')

    for line in dnaFileIN:
        dnaSeq = line[:-1]

        for i in range(len(dnaSeq)):

            # Generate key to categorie other sequences which call into this template
            dnaKey = dnaSeq[:i] + '_' + dnaSeq[i+1:]

            if dnaKey in dnaDict:
                # Add the sequence as a value to this key
                dnaDict[dnaKey].append(dnaSeq)
            # Add as a new key
            else:
                dnaDict[dnaKey] = [dnaSeq]
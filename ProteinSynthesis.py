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

# Traverse the cell grouping board
# Find moves to allow knight to visit every square on the board once
def cellGroupGraph(bdSize):

    cgGraph = AL_Graph()

    for row in range(bdSize):
        for col in range(bdSize):
            # Get a list index for the (row, col) value
            currNodeIndex = convertCoord(row, col, bdSize)
            # Get all positions this knight could potentially move to
            newPositions = genLegalPositions(row, col, bdSize)

            for position in newPositions:
                newNodeIndex = convertCoord(position[0], position[1], bdSize)
                cgGraph.addEdge(currNodeIndex, newNodeIndex)

    return cgGraph

# Convert a x,y coordinate into an index for list
def convertCoord(x, y, bdSize):

    index = (x * bdSize) + y

    return index

def genLegalPositions(x, y, bdSize):

    newPositions = []

    moveRules = [(-1, -2), (-1, 2), (-2, -1), (-2,  1),
                 ( 1, -2), ( 1, 2), ( 2, -1), ( 2, -1)]

    # Generate all possible new positions from th current position
    for move in moveRules:

        newX = x + move[0]
        newY = y + move[1]

        # Check if new position is actually possible on this board
        if legalCoord(newX, bdSize) and legalCoord(newY, bdSize):

            newPositions.append((newX, newY))

    # Return all positions we could potentially move to
    return newPositions

def legalCoord(coord, bdSize):

    if coord >= 0 and coord <= bdSize:
        return True

    else:
        return False
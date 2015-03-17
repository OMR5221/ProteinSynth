import sys


#In this problem, we ask a simple question:

#How many times can one string occur as a substring of another?
#Recall from Frequent Words Problem that different occurrences of a substring can overlap with each other.

#For example, ATA occurs three times in CGATATATCCATAG

#Pattern Matching Problem

#Find all occurrences of a pattern in a string.

#Given: Strings Pattern and Genome.

#Return: All starting positions in Genome where Pattern appears as a substring. Use 0-based indexing.

#Sample Dataset (Input)

# Pattern: ATAT
#Genome: GATATATGCATATACTT

# Sample Output: Start Index in Genome

# 1 3 9

patternStr = "ATAT"
genomeStr = "GATATATGCATATACTT"

pattern_StartPos = []

# Using the Python function


# Looping through String
stopPatternSearch = False
stopGenomeSearch = False

patternList = list(patternStr)
genomeList = list(genomeStr)

startGenomeIndex = 0
currGenomeIndex = 0
patternIndex = 0

patternComplete = False


while startGenomeIndex < len(genomeList):

    #print("Start Genome Index: " + str(startGenomeIndex) + "\n")

    patternIndex = 0
    stopPatternSearch = False
    patternComplete = False

    currGenomeIndex = startGenomeIndex

    while patternIndex < len(patternList) and not stopPatternSearch:

        #print("\tStart Pattern Index: " + str(startGenomeIndex) + "\n")

        if genomeList[currGenomeIndex] == patternList[patternIndex]:

            if patternIndex == len(patternStr)-1:
                patternComplete = True
                stopPatternSearch = True
                currGenomeIndex += 1

            else:
                currGenomeIndex += 1
                patternIndex += 1
        else:
            stopPatternSearch = True

    if patternComplete:
        pattern_StartPos.append(startGenomeIndex)

    startGenomeIndex += 1

print(pattern_StartPos)
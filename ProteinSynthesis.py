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

class BioServer:

    def __init__(self):
        serverName = ""
        driver = "MySQLdb"
        user = "root"
        passwd = "gtrhalo2"
        host = "127.0.0.1"
        db = "bioseqdb"
        connected = False

    def get_Connected(self):
        return self.connected

    def set_Connnected(self, connectAction):
        self.connected = connectAction

    def Connect_Server(self, inDriver, inUser, inPasswd, inHost, inDB):

        if self.connected == False:
            print("Should the default settings be used to connect? [Y/N]")
            defaultConnInput = raw_input()

            if defaultConnInput == "Y" or defaultConnInput == "y":
                server = BioSeqDatabase.open_database(driver=self.driver, user=self.user, passwd=self.passwd, host=self.host, db=self.db)
            else:
                print("Please enter the following BioSQL Server settings: \n")
                print("\tDriver: \n")
                inDriver = raw_input()
                print("\tUsername: \n")
                inUser = raw_input()
                print("\tPassword: \n")
                inPasswd = raw_input()
                print("\tHost Name: \n")
                inHost = raw_input()
                print("\tDatabase Name: \n")
                inDB = raw_input()
                print("Attempting Connection: \n" )
                server = BioSeqDatabase.open_database(driver=self.driver, user=self.user, passwd=self.passwd, host=self.host, db=self.db)

        else:
            print("Already connected to a server!")

    def Create_SubDB(self):
    # 1. Create a sub database for sequences to be added to
    server = BioSeqDatabase.open_database(driver="MySQLdb", user="root", passwd="gtrhalo2", host="127.0.0.1", db="bioseqdb")

    db = server.new_database("orchids", description="Just to Test")

    def Connect_DB(self):

        # Fetch the wanted sequences and load into the created sub database
        db = server["orchids"]

    def Connect_SvrDB(self):
            server = BioSeqDatabase.open_database(driver="MySQLdb", user="root", passwd="gtrhalo2", host="127.0.0.1", db="bioseqdb")

            # Fetch the wanted sequences and load into the created sub database
            db = server["orchids"]

    def FetchSeq(self):

        if connectedToServer:
            # Get server Already connected to
        else:
            server = Connect_Server()

        if connectedToDB:
            # Get db connected to
        else:
            db = Connect_DB()

        # Fetch the following records from the Web
        handle = Entrez.efetch(db="nuccore", id="6273291,6273290,6273289", rettype="gb", retmode="text")

        # Load returns number of records handled into db
        count = db.load(SeqIO.parse(handle, "genbank"))

        print "Loaded %i records" % count

        handle.close()

        # Commit records to db to ensure they are entered
        server.adaptor.commit()

    def FetchSeq(self):
        # Extract sequences from sub database for use
        # Connect to "server" and open db
        server = BioSeqDatabase.open_database(driver="MySQLdb", user="root", passwd="gtrhalo2", host="127.0.0.1", db="bioseqdb")

        # select sub database
        db = server["orchids"]

        print("This database contains %i records" % len(db))

        # find records for each identifier
        for identifier in ['6273290']:
            # Retrieve the sequence record by lookuping the Id
            seq_record = db.lookup(gi=identifier)

            print(seq_record.id, seq_record.description[:50] + "...")

            print("Sequence Length: %i\n" % len(seq_record.seq))

    def Disconnect(self):
        self.connected = False

#Structure of each cells functions, time and quantity
class Cell:

    def __init__(self, id, dnaSeq, maxPSNum, maxMutNum):
        self.id = id
        self.startDNASeq = dnaSeq
        self.currDNASeq = MutableSeq(dnaSeq, IUPAC.unambiguous_dna)
        self.mrnaSeq
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
        self.mrnaSeq = coding_strand.transcribe()

        # TRue Biological Process
        mRNA_strandB = template_strand.reverse_complement().transcribe()


    def Translation(self):
        # Find Proteins attributed to a mrnaSeq
        messenger_rna = Seq(self.mrnaSeq, IUPAC.unambiguous_rna)

        # Stop translation at first in frame stop codon as described in dna table
        messenger_rna.translate(to_stop=True)

        # We can also specify our own stop codon
        # Use for our Simulation by splicing dna seq at various points?
        messenger_rna.translate(stop_symbol="@")


    def ProteinSynthesis(self):

        self.rnaSeq = []
        self.numPS += 1

        # Replication of dna strand into sense (coding) strand and anti-sense strand (template)

        # Transcribe dna strand and find complement rna sequence via RNA polymerase

    # Define Death of a cell by (changing graph color and removing resources and setting weight to 0)
    def Apoptosis(self):
        None

    # Mutate dna sequence by lering nucleotides randomly
    def DNA_Damage(self, mutateNum):
        '''
        if self.mt_Num == mutateNum:
            # Alter some nucleotides
        '''
        None

    def DNA_Repair(self):
        None

    def __iter__(self):
        return self.__iter__()

def ProteinSynthesis_Sim(numSeconds, dnaSeq, maxPSNum, numCells):

    ps_cellLists = []

    for i in range(maxPSNum):

        ps_cellLists.insert(i, [])

    for i in range(numCells):

        aCell = Cell(i, dnaSeq, maxPSNum)
        '''
        if ps_cellLists[aCell.ps_srtNum] == None:

            ps_cellLists[aCell.ps_srtNum] = []
            ps_cellLists[aCell.ps_srtNum].append(aCell)

        else:
        '''
        ps_cellLists[aCell.ps_srtNum].append(aCell)

    ''' Create a balance tree of cells organized by process_StartNum
    rootNum = maxPSNum / 2
    rootNode = Tree(rootNum)

    psTree = BinarySearchTree()
    psTree.root = rootNode
    '''

    for currentSecond in range(numSeconds):

        #  Check if a cell is ready to perform protein synthesis
        genTime = GenPSTime(maxPSNum)

        # Run all cells whose startTime is less than the time Number generated
        for i in range(genTime):
            for cell in ps_cellLists[i]:
                cell.ProteinSynthesis()
                # print("Cell " + str(cell.id) + ": " + str(cell.rnaSeq))

    for i in range(maxPSNum):
        for cell in ps_cellLists[i]:
            print("Cell " + str(cell.id) + "(" + str(cell.numPS) + "): " + str(cell.rnaSeq))

def GenPSTime(maxPSNum):

    return random.randrange(1, maxPSNum+1)



# Create A Menu to display
showMenu = True
menuSelection = 0

myBioServer = BioServer()

myCell = Cell()

while showMenu:
    print("Welcome to Protein Synthesis Simulation!\n")
    print("Please select an option below.\n\n")
    print("1. Connect to a BioSQL Server\n")
    print("2. Create a Sub Database\n")
    print("3. Connect to a Sub Database\n")
    print("4. Fetch and Load a Sequence\n")
    print("5. Run Protein Synthesis\n")

    raw_input(menuSelection)

    if menuSelection == 1:
        if myBioServer.connected == False:
            print("Initiating BioSQL Server Connection Process:\n")
            myBioServer.Connect()
        else:
            print("You are already connected to a BioSQL Server: \n")
            print("\t%s as user %s" % (myBioServer.serverName, myBioServer.user))

    elif menuSelection == 2:
        myBioServer.Create_SubDB()

    elif menuSelection == 3:
        myBioServer.Connect_Server()
        myBioServer.Connect_DB()

    elif menuSelection == 4:
        myBioServer.Connect_Server()
        myBioServer.Connect_DB()
        myCell.currDNASeq = myBioServer.FetchSeq()

    elif menuSelection == 5:
        ProteinSynthesis_Sim(100, str(seq_record.seq), 200, 15)

    else:
        print("Error")


# Create cells and run program to simulate
# program has pieces which move in various ways through the knights tour problem
# but are dictated by AI logic as to which path to take

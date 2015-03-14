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
        driverName = "MySQLdb"
        userName = "root"
        passwdName = "gtrhalo2"
        hostName = "127.0.0.1"
        dbName = "bioseqdb"
        connected = False
        server = None
        db = None

    def get_ConnectStatus(self):
        return self.connected

    def set_ConnnectStatus(self, connectAction):
        self.connected = connectAction

    def Connect_Server(self, inDriver, inUser, inPasswd, inHost, inDB):

        if self.connected == False:
            print("Should the default settings be used to connect? [Y/N]")
            defaultConnInput = raw_input()

            if defaultConnInput == "Y" or defaultConnInput == "y":
                self.server = BioSeqDatabase.open_database(driver=self.driverName, user=self.userName, passwd=self.passwdName, host=self.hostName, db=self.dbName)
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
                self.server = BioSeqDatabase.open_database(driver=inDriver, user=inUser, passwd=inPasswd, host=inHost, db=inDB)

        else:
            print("Already connected to a server!")

    def Create_SubDB(self):

        self.db = self.server.new_database("orchids", description="Just to Test")

    def Connect_DB(self):

        # Fetch the wanted sequences and load into the created sub database
        self.db = self.server["orchids"]

    def FetchSeq(self):

        if not self.connected:
            self.Connect_Server()

        if self.db == None:
            self.Connect_DB()

        # Fetch the following records from the Web
        handle = Entrez.efetch(db="nuccore", id="6273291,6273290,6273289", rettype="gb", retmode="text")

        # Load returns number of records handled into db
        count = self.db.load(SeqIO.parse(handle, "genbank"))

        print "Loaded %i records" % count

        handle.close()

        # Commit records to db to ensure they are entered
        self.server.adaptor.commit()

    def FetchSeq(self):
        # Extract sequences from sub database for use
        # Connect to "server" and open db
        # server = BioSeqDatabase.open_database(driver="MySQLdb", user="root", passwd="gtrhalo2", host="127.0.0.1", db="bioseqdb")

        # select sub database
        # db = server["orchids"]

        print("This database contains %i records" % len(self.db))

        # find records for each identifier
        for identifier in ['6273290']:
            # Retrieve the sequence record by lookuping the Id
            seq_record = self.db.lookup(gi=identifier)

            print(seq_record.id, seq_record.description[:50] + "...")

            print("Sequence Length: %i\n" % len(seq_record.seq))

    def Disconnect(self):
        self.connected = False

# Create A Menu to display
showMenu = True
menuSelection = 0

myBioServer = BioServer()

while showMenu:
    print("Welcome to Protein Synthesis Simulation!\n")
    print("Please select an option below.\n\n")
    print("1. Connect to a BioSQL Server\n")
    print("2. Create a Sub Database\n")
    print("3. Connect to a Sub Database\n")
    print("4. Fetch and Load a Sequence\n")
    print("5. Run Protein Synthesis\n")

    biosql_menuSelect = int(raw_input())

    if biosql_menuSelect == 1:

        print("\t1. Connect to a BioSQL Server\n")
        print("\t2. Create a Sub Database\n")
        print("\t3. Connect to a Sub Database\n")
        print("\t4. Fetch and Load a Sequence\n")

        if biosql_menuSelect == 1:

            if myBioServer.connected == False:
                print("Initiating BioSQL Server Connection Process:\n")
                myBioServer.Connect()
            else:
                print("You are already connected to a BioSQL Server: \n")
                print("\t%s as user %s" % (myBioServer.serverName, myBioServer.userName))

        elif biosql_menuSelect == 2:
            myBioServer.Create_SubDB()

        elif biosql_menuSelect == 3:
            myBioServer.Connect_Server()
            myBioServer.Connect_DB()

        elif biosql_menuSelect == 4:
            myBioServer.Connect_Server()
            myBioServer.Connect_DB()
            myCell.currDNASeq = myBioServer.FetchSeq()

        elif biosql_menuSelect == 5:
            ProteinSynthesis_Sim(100, str(seq_record.seq), 200, 15)

    else:
        print("Error")

# Create cells and run program to simulate
# program has pieces which move in various ways through the knights tour problem
# but are dictated by AI logic as to which path to take

for seqRecord in SeqIO.parse("ls_orchid.fasta", "fasta"):
    # print(seqRecord.id)
    #print(seqRecord.seq)
    ProteinSynthesis_Sim(500, str(seqRecord.seq), 200, 5, 200)
    # print(len(seqRecord))

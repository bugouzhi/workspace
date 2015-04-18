#!/usr/bin/python
#a script to extract all identified proteins from a search result file
import sys
import getopt
from sets import Set

def parseSearchResult(searchResult, protInd):
    proteins = Set();
    results = open(searchResult, 'r');
    for result in results:
        result.strip();
        tokens = result.split("\t");
        proteins.add(tokens[protInd]);
    #print "IDed total proteins " + str(len(proteins));
    return proteins;

def parseFastaFile(fasta):
    proteinMap = {};
    fastaFile = open(fasta, 'r');
    protein = "";
    proteinName = "";
    for line in fastaFile:
        line = line[:-1];
        if line.startswith('>'):
            if proteinName is not "":
                proteinMap[proteinName] = protein;
                protein = "";
            proteinName = line[1:];
        else:
            protein = protein + line;
    return proteinMap;

def getIdedProteins(searchResult, fastaFile, index, outfile):
    out = open(outfile, 'w');
    identified = parseSearchResult(searchResult, index);
    print "identified " + str(len(identified));
    proteins = parseFastaFile(fastaFile);
    print "total: " + str(len(proteins));
    for ID in identified:
        for prot in proteins:
            if(prot.startswith(ID)):
                out.write(">" + prot + "\n" + proteins[prot] + "\n");
    
def usage():
    print "<search result> <fasta file> <protein column index> <outfile>"    

def main():
    usage();
    getIdedProteins(sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]);
    


if __name__ == "__main__":
    main()

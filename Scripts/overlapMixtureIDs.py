#!/usr/bin/python
import sys
import getopt
from sets import Set




def parseResultFile(resultFile):
    results = open(resultFile, 'r');
    resultList = [];
    for line in results:
        #line.rstrip();
        line = line[:-1];
        tokens = line.split(" & ");
        resultSet = set(tokens);
        resultList.append(resultSet);
    return resultList;    

def usage():
    print "<mixIDs1> <mixIDs2>" 
        
def main():
    usage();
    resultList1 = parseResultFile(sys.argv[1]);
    resultList2 = parseResultFile(sys.argv[2]);
    for index in range(len(resultList1)):
        result1 = resultList1[index];
        maxoverlap=0;
        for index2 in range(len(resultList2)):
            result2 = resultList2[index2];
            #print "set 1 " + str(result1);
            #print "set 2 " + str(result2);
            overlap = len(result1 & result2);
            if(overlap > maxoverlap):
                maxoverlap = overlap;
        print "" + str(index) + "\toverlapIDs:\t" + str(maxoverlap); 


if __name__ == "__main__":
    main()


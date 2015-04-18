#!/usr/bin/python
import sys
import glob
import subprocess
import os

def main():
    result_file = sys.argv[1];
    out_file = sys.argv[2];
    cmd = "cat " + result_file + " | awk -F '\t' '{if($15 < 0.01) print}' | sort -t '\t' -k8,8 -k7,7n -k12,12g | sort -t '\t' -k8,8 -k7,7n -u > " + out_file;
    path = os.getcwd();
    #cmd = "more " + path + "/" + result_file + "\n";
    #cmd = "cat " + result_file;
    print "command " + cmd;
    #print "current path " + os.getcwd();
    #os.system(cmd);
    subprocess.call(cmd, shell=True);
    print "current path " + path;
    #subprocess.Popen(cmd, stdout=PIPE).stdout.read();
    #subprocess.call(cmd, shell=True);
    #subprocess.check_output(cmd);

if __name__ == "__main__":
    main()

#!/usr/bin/env python
import os, re
import commands
import math, time
import sys
import optparse
import ROOT
import pickle
print 
print 



def main():

    #configuration                                                                                                               
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--in',  dest='input', help='input directory', default=None,  type='string')
    (opt, args) = parser.parse_args()

    queue = "1nh"
    from subprocess import Popen, PIPE
    directory=opt.input
    print 'looking into: '+directory+'...'
    prepend='root://eoscms.cern.ch//eos/cms'

    path = os.getcwd()
    print
    print 'do not worry about folder creation:'
    os.system("rm -r tmp")
    os.system("mkdir tmp")
    os.system("mkdir res_test")
    print



    eos_cmd = '/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select'
    data = Popen([eos_cmd, 'ls', '/eos/cms/'+directory],stdout=PIPE)
    out,err = data.communicate()
    full_list = []
    for line in out.split('\n'):
        if len(line.split()) == 0: continue
        full_list.append(prepend + directory + '/' + line)


    input_list=[]
    input_list=full_list
    x=1
    for ifile in xrange(0,len(input_list)):
        inF=input_list[ifile]
        outF="output_"+str(x)
        os.system("mkdir tmp/"+str(x))
        os.chdir("tmp/"+str(x))
        print "input file : " , inF
        with open('job.sh', 'w') as fout:
            fout.write("#!/bin/sh\n")
            fout.write("echo\n")
            fout.write("echo\n")
            fout.write("echo 'START---------------'\n")
            fout.write("echo 'WORKDIR ' ${PWD}\n")
            fout.write("cd "+str(path)+"\n")
            fout.write("root -l -b <<EOF\n")
            fout.write("EOF \n")
            fout.write("gSystem->Load("'"libFWCoreFWLite.so"'") \n")
            fout.write("AutoLibraryLoader::enable() \n")
            fout.write("gSystem->Load("'"libUserCodeTopAnalysis.so"'") \n")
            fout.write(".L ReadTree.cc+  \n")
#            fout.write("ReadTree(" '"'+inF+'"' ",res_"'"'+str(x)+'"'") \n")
#            fout.write("ReadTree(" '"'+inF+'"' '"'+outF +'"' " ) \n")
            fout.write("ReadTree(" '"'+inF+'"' "," '"'+outF +'"' " ) \n")
            fout.write(".q; \n")
            fout.write("EOF \n")
            fout.write("EOF \n")
            fout.write("echo 'STOP---------------'\n")
            fout.write("echo\n")
            fout.write("echo\n")
        os.system("chmod 755 job.sh")
   
   ###### sends bjobs ######
#        os.system("sh job.sh")
        os.system("bsub -q "+queue+" -o logs job.sh")
        #print "job nr " + str(x) + " submitted"
        os.chdir("../..")
        print
        print "your jobs:"
        print "output file : ", outF
        #    os.system("bjobs")
        print
        print 'END'
        print
        x = x+1
        print " value of x", x
if __name__ == "__main__":
    sys.exit(main())

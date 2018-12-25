

#!/usr/bin/python
 
import sys
import re
 
class Compatibility:
  pass
args = Compatibility
 
if sys.version_info < (2,7):
  if len(sys.argv)<2:
    print("splitLHE.py requires 3 command line arguments, inputFile ouFileNameBase and nFiles" )
    sys.exit(1)
  args.inputFile = sys.argv[1]
#  args.outFileNameBase = sys.argv[2]
#  args.nFiles = int(sys.argv[3])
else:
  import argparse
  parser = argparse.ArgumentPaerser(descrition="Splits LHE events from input file into N output files")
  parser.add_argument("inputFile",help="Input LHE file name")
  parser.add_argument("outFileNameBase",help="Output LHE file name Base (a number and .he will be added to it)")
  parser.add_argument("nFiles",help="Number of files to split the events between",type=int)
 
#  args = parser.parse_args()
 
#if args.nFiles < 2:
#  print("Error: nFiles must be > 1")
#  sys.exit(1)
 
fin = ""
try:
  fin = open(args.inputFile)
except:
  print("Error: Input file: %s could not be opened, exiting." % args.inputFile)
  sys.exit(1)
 
eventNum = 0
init = False
inFooter = False
footLines = []
for line in fin:
  if re.match(r"[^#]*</LesHouchesEvents>",line):
    inFooter = True
    footLines.append(line+"\n")
  elif inFooter:
    footLines.append(line+"\n")
  elif init:  
    if re.match(r"[^#]*</event>",line):
      eventNum += 1
  elif re.match(r"[^#]*</init>",line):
    init = True
 
eventsTotal = eventNum
print "N Events Total: %i" % eventsTotal



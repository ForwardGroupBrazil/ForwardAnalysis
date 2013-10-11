#!/usr/bin/python

import sys, getopt, os

def main(argv):

   user_ = "dmf"
   try:
      opts, args = getopt.getopt(argv,"hu:",["user="])
   except getopt.GetoptError:
      print 'CleanCrabServers.py -u <number of jobs>'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print 'CleanCrabServers.py -u <number of jobs>'
         sys.exit()
      elif opt in ("-u", "--user"):
         user_ = arg

   print ""
   print "=========="
   print "Clean Tool"
   print "=========="
   print ""
   listOfservers = ['hcc-crabserver.unl.edu','vocms20.cern.ch','vocms83.cern.ch','submit-4.t2.ucsd.edu','submit-6.t2.ucsd.edu']

   for line in listOfservers:
     os.environ["SERVERNAME"] = line
     print "<<< SERVER ", line, " >>>"
     print "- List of Files: "
     os.system("gsissh $SERVERNAME \"ls -lgort\"")
     command = "gsissh $SERVERNAME \"rm -rf "+ user_ + "*\""
     #print command
     #readuser = os.environ["USER"]
     #print readuser
     os.system(command)

   print ""
   print "===="
   print "End!"
   print "===="
   print ""

if __name__ == "__main__":
   main(sys.argv[1:])


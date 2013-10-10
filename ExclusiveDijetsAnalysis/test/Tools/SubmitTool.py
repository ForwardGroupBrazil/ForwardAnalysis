#!/usr/bin/python

import sys, getopt, os

def main(argv):
   numj = 5010
   try:
      opts, args = getopt.getopt(argv,"hn:",["njobs="])
   except getopt.GetoptError:
      print 'SubmitTool.py -n <number of jobs>'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print 'SubmitTool.py -n <number of jobs>'
         sys.exit()
      elif opt in ("-n", "--njobs"):
         numj = arg

   print ""
   print "==========="
   print "Submit Tool"
   print "==========="
   print ""

   total=int(numj)/500

   count = 0
   while count < total+1:
      os.system("multicrab -submit 500")
      count=count+1
      print "\nInteraction",count,"from total",total+1
   else:
      print "\nAll tasks submitted."

   print ""
   print "===="
   print "End!"
   print "===="
   print ""

if __name__ == "__main__":
   main(sys.argv[1:])

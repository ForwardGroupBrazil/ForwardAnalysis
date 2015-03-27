#!/usr/bin/python

import sys, getopt, os

def main(argv):
   numj = 500
   pathname = ''
   try:
      opts, args = getopt.getopt(argv,"hn:c:",["njobs=","path="])
   except getopt.GetoptError:
      print 'SubmitTool.py -n <number of jobs> -c <path>'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print 'SubmitTool.py -n <number of jobs> -c <path>'
         sys.exit()
      elif opt in ("-n", "--njobs"):
         numj = arg
      elif opt in ("-c", "--path"):
         pathname = arg
   print ""
   print "==========="
   print "Submit Tool"
   print "==========="
   print ""

   total=int(numj)/500
   count = 0
   
   if(len(pathname) > 0):
     command = "multicrab -submit 500 -c "+pathname 
   else:
     command = "multicrab -submit 500"
   while count < total+1:
      os.system(command)
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

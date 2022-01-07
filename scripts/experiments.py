# Script to set up parameters for MibS instance path.
# Used to run experiments for diffrent cuts 
# Last edited by yux616
# Apr 2020

# Executable path and name
pbsfile = '/home/ted/Projects/Cbc/scripts/batch.pbs'

# Instance path
# Directory name and path containing test instances in .mps format
# Keys are used to name subdirs in output dir
testSets = {
    'MIPLIB2017': '/home/ted/DataSets/miplib2017'
}

versions = {
#    'before' : '/home/ted/Projects/build-cbc-before/bin/cbc',
    'after' : '/home/ted/Projects/build-cbc-master/bin/cbc'
}

# Output parent path
outputDir = '/home/ted/Projects/Cbc/output'

# Name
testname = 'cbc'

# Set up senarios
scenarios = {}

commonParams = {
    'sec': '3600',
    '-solve' : ''
}    

scenarios['default'] = {

}
    

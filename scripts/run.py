# Script to run IBLP
# The script may produce auxiliary folder/files in run directory. 
# Some os function requires Python 3.5+

import os
import shutil, subprocess
import matplotlib.pyplot as plt
from experiments import testSets, outputDir, scenarios, pbsfile, testname, versions, commonParams

def runExperiments(testSets, outDir, versions, scenarios):
    """
        Use to run experiments on local machine. 
    """

    # set up output directories
    # use hierarchy:  outDir/version/scenario_name/testset_name/
    for v in versions:
        for s in scenarios:
            currpath = os.path.join(outDir, v, s)
            if not os.path.exists(currpath):
                os.mkdir(currpath)
            for t in testSets:
                currsubpath1 = os.path.join(currpath, t)
                if not os.path.exists(currsubpath1):
                    os.mkdir(currsubpath1)
   
    # run experiments use command line paramters
    for v in versions:
        exe = versions[v]
        for s in scenarios:
            for t in testSets:
                outsubpath = os.path.join(outDir, v, s, t)
                os.chdir(outsubpath)   
                with os.scandir(testSets[t]) as inst_it: 
                    for instance in inst_it:
                        if not (instance.name.endswith('.mps')
                                or instance.name.endswith('.mps.gz')):
                            continue
                        for p in scenarios[s]:
                            params.append(p)
                            params.append(scenarios[s][p])
                        if instance.name.endswith('.mps'):
                            outname = instance.name[:-4]+'.out'
                        else:
                            outname = instance.name[:-7]+'.out'
                        if os.path.exists(outname):
                            print("File", outname, "exists!")
                            continue
                        print ("Now processing ", instance.name)
                        print ("Args: ", params[2:])
                        outfile = open(outname,'w')
                        output = subprocess.check_output(params, universal_newlines = True)
                        print (output)
                        outfile.write(output)
                        outfile.close()
    return

def runExperimentsPBS(testSets, outDir, versions, scenarios, pbsfile):
    """
        Use to submit batch jobs via qsub. 
    """
    # set up output directories
    # use hierarchy:  outDir/version/param_scenario_name/testset_name
    for v in versions:
        for s in scenarios:
            for t in testSets:
                currsubpath = os.path.join(outDir, v, s, t)
                if not os.path.exists(currsubpath):
                    os.makedirs(currsubpath)
                else:
                    print('Directory already exists! Ignoring...')

    # submit experiments use command line paramters
    # see .pbs file anc comments for job submission arguments
    for v in versions:
        exe = versions[v]
        for s in scenarios:
            for t in testSets:
                outsubpath = os.path.join(outDir, v, s, t)
                os.chdir(outsubpath)   
                with os.scandir(testSets[t]) as inst_it: 
                    for instance in inst_it:
                        if not (instance.name.endswith('.mps')
                                or instance.name.endswith('.mps.gz')):
                            continue
                        params = ''
                        for p in scenarios[s]:
                            params +=  (p + ' ' + scenarios[s][p] + ' ')
                        outfile = os.path.join(outsubpath, instance.name[:-4]+'.out')
                        errfile = os.path.join(outsubpath, instance.name[:-4]+'.err')
                        if os.path.exists(outfile):
                            print("File", outfile, "exists!")
                            continue
                        subprocess.run(["qsub", "-v", 
                                        "EXECUTABLE="+exe+","
                                        +"INSTANCENAME="+instance.path+","
                                        +"PARAMARG="+params,
                                        "-o", outfile,
                                        "-e", errfile,
                                        "-N", testname,
                                        pbsfile])
    return                    

if __name__ == "__main__":

    for s in scenarios:
        scenarios[s].update(commonParams)

    ######################### Run Experimests #########################
    # local: provide paths in runparams.py
    # runExperiments(testSets, outputDir, versions, scenarios)
    
    # using pbs file: provide paths in runparams.py
    runExperimentsPBS(testSets, outputDir, versions, scenarios, pbsfile)
    



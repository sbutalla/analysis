#############################################
## Script to get production cross-section, ##
## branching ratio, and kinetmatic and     ##
## geometric acceptances. Interfaces with  ##
## MadGraph5 and MadAnalysis5.             ##
## Stephen D. Butalla                      ##
## 2022/02/13                              ##
#############################################

import argparse
import gzip
import numpy as np
import os
import sys
from time import sleep
from tqdm import tqdm

# particle mass parameters
muonMass        = 0.1 # GeV
default_fd1Mass = 5   # GeV
default_fd2Mass = 2   # GeV

# relative directory paths and file names
baseMG5Dir     = "/afs/cern.ch/user/s/sbutalla/analysis/Madgraph/MG5_aMC_v2_6_7"
massScanDir    = "massScans" # results dir for this script
modelDir       = "models/DMsimp" # UFO file dir
outputDir      = "DMsimp" # output dir after process generated
outputEventDir = "DMsimp/Events/run_01" # lhe file dir 
ma5Dir         = "ma5_test/madanalysis5" # ma5 base dir
ma5_lheDir     = "ma5_test/madanalysis5/DMsimp" # ma5 lhe file dir
ma5_script     = "cuts_script.txt" # ma5 analysis script
ma5_latex      = "Output/DVI/MadAnalysis5job_0" # ma5 LaTeX output dir

# command extensions for quiet mode
suppress = " > /dev/null 2>&1" # forwards stdout and stderr to /dev/null 

# production process and fill process
processes = ["p p > Zd", "p p > Zd, (Zd > FD11 FD11, (FD11 > FD21 mu+ mu-))"]

class colors:
    MAGENTA   = '\033[95m'
    BLUE      = '\033[94m'
    CYAN      = '\033[96m'
    GREEN     = '\033[92m'
    YELLOW    = '\033[93m'
    RED       = '\033[91m'
    ENDC      = '\033[0m'
    BOLD      = '\033[1m'
    UNDERLINE = '\033[4m'

def newScanSetup(runName):
    #####################################################
    # Function that creates the new mass scan directory #
    # if it is the first time running this software,    #
    # and creates the run subdirectory.                 #
    # Arguments:                                        #
    # runName     - String; run name set by the user    #
    # Returns:                                          #
    # fullPath - String; relative path to the results   #
    #            directory                              #
    #####################################################
    try:
        os.makedirs(massScanDir) # create directory for mass scan runs
    except FileExistsError: # skip if directory already exists
        pass
    
    fullPath = massScanDir + "/" + runName
    try:
        os.makedirs(fullPath) # create directory for specific run
    except FileExistsError: # skip if directory already exists
        pass
    
    return fullPath

def getRecentDir():
    #####################################################
    # Function that retrieves the most recently         #
    # created directory.                                #
    # Arguments:                                        #
    # dirPath - String; path to the base dir where sub  #
    #           dirs are located                        #
    # Returns:                                          #
    # targetDir - String; most recently created dir     # 
    #####################################################
    print("Getting most recent dir at " + baseMG5Dir + "/" + ma5Dir)
    baseDir = baseMG5Dir + "/" + ma5Dir
    folders = next(os.walk(baseDir))[1]
    
    time = np.empty(len(folders))
    counter = 0
    for folder in range(len(folders)):
        time[counter] = os.path.getctime(baseDir + "/" + folders[folder])
        counter += 1

    recentDir = np.argmax(time)

    return folders[recentDir]

def prepareDataFiles(fullPath, prodCrossSection, branchingRatio, acceptances):
    #####################################################
    # Function that initializes the data files (.txt)   #
    # format. Returns nothing.                          #
    # Arguments:                                        #
    # fullPath         - String; (relative) path to the #
    #                     base dir where sub dirs are   #
    #                     located                       #
    # prodCrossSection - Boolean; compute prod. xsec    #
    # branchingRatio   - Boolean; compute BR            #                                  
    # acceptances      - Boolean; compute geo. + kin.   #
    #                    acceptances                    #
    #####################################################
    if prodCrossSection == True:
        fileName = fullPath + "/prod_crossX.txt"
        fileCrossSection = open(fileName,"w+")
        fileCrossSection.write("ZD_mass    numEvents    prod_crossX\n")
        fileCrossSection.close()
    else:
        pass
    
    if branchingRatio == True:
        fileName = fullPath + "/branching_ratio.txt"
        fileBR = open(fileName,"w+")
        fileBR.write("ZD_mass    fd1_mass    fd2_mass    ZD_tofd1fd1bar_BR    fd1_tofd2_diMuons_BR    processBR\n")
        fileBR.close()
    else:
        pass

    if acceptances == True:
        fileName = fullPath + "/acceptances.txt"
        fileAcc = open(fileName,"w+")
        fileAcc.write("ZD_mass    fd1_mass    fd2_mass    cumulative_acceptance\n")
        fileAcc.close()
    else:
        pass

def editParamCard(ZDmass, fd1Mass, fd2Mass, modelDir):
    #####################################################
    # Function that edits the parameters.py file.       #
    # Returns nothing.                                  # 
    # Arguments:                                        #
    # ZDmass   - Float; mass of dark photon (ZD)        #
    # fd1Mass  - Float; mass of dark fermion 1 (fd1)    #
    # fd2Mass  - Float; mass of dark fermion 2 (fd2)    #                                  
    # modelDir - String; directory where model info is  #
    #            (FeynRules UFO format)                 #
    #####################################################
    fileName = modelDir + "/parameters.py"
    with open(fileName, 'r') as file: # open and read file
        contents = file.readlines()

    contents[254] = "                value = %d,\n" % fd1Mass # Replace fd1 mass in the parameters.py file
    contents[262] = "                value = %d,\n" % fd2Mass # Replace fd1 mass in the parameters.py file
    contents[270] = "                value = %d,\n" % ZDmass  # Replace ZD mass in the parameters.py file

    with open(fileName, 'w') as file: # overwrite file
        file.writelines(contents)

def editProcCard(process, computeWidths = False):
    #####################################################
    # Function that edits the MG5 process card          # 
    # (proc_card.dat). Returns nothing.                 #
    # Arguments:                                        #
    # process       - String; process for MG5 to        #
    #                 generate (must be in MG5 format)  #
    # computeWidths - Boolean; Select whether to        #
    #                 compute the decay widths of the   #
    #                 process. Default = False          #
    #####################################################
    fileName = "proc_card.dat"
    with open(fileName, 'r') as file:
        contents = file.readlines()
    
    contents[25] = "generate %s" % process + "\n" # Replace process in proc_card

    if computeWidths == True:
        contents[26] = "compute_widths all --body_decay=4\n"
    elif computeWidths == False:
        contents[26] = "\n"

    with open(fileName, 'w') as file: # overwrite file
        file.writelines(contents)

def getXsec(quiet):
    #####################################################
    # Function that retrieves the num. of events prod.  #
    # xsec. from the lhe file after event generation.   #
    # Returns:                                          #
    # numEvents - Int; the number of simulated events   #
    # prodXsec  - Float; the prod. xsec                 #
    # Values returned as seperate variables.            #
    #####################################################
    fileName = outputEventDir + "/unweighted_events.lhe.gz"

    if quiet == True:
        os.system('gzip -d ' + fileName + suppress) # gunzips lhe file; can't parse the file without decompressing
    else:
        os.system('gzip -d ' + fileName)
    
    decompressedFile = fileName[:len(fileName)-3] # get file name of decompressed lhe file
    
    with open(decompressedFile, 'r') as file: # open and read file
        contents = file.readlines()
    
    numEvents = int(contents[523].split()[-1]) # get the number of events
    prodXsec = float(contents[524].split()[-1]) # get the production cross-section
    #numEvents = int(contents[543].split()[-1]) # get the number of events
    #prodXsec = float(contents[544].split()[-1]) # get the production cross-section
    
    return numEvents, prodXsec

def getBR():
    #####################################################
    # Function that retrieves the partial and total BRs #
    # from the parameter card (param_card.dat). Accepts #
    # no args.                                          #
    # Returns:                                          #
    # ZD_tofd1fd1bar_BR    - Float; partial BR of       #
    #                        ZD -> fd1 fd1~             #
    # fd1_tofd2_diMuons_BR - Float; partial BR of       #
    #                        fd1 -> fd2 mu+ mu-         #
    # processBR            - Float; total process BR    #
    # Values returned as seperate variables.            #
    #####################################################
    fileName = modelDir + "/param_card.dat"
    with open(fileName, 'r') as file:
        contents = file.readlines()
        
    ZD_tofd1fd1bar_BR = float(contents[174].split()[0]) # get BR for ZD to fd1/fd1~
    fd1_tofd2_diMuons_BR = float(contents[192].split()[0]) # get BR for fd1 ro fd2/fd2~ and dimuons
    
    processBR = ZD_tofd1fd1bar_BR * np.power(fd1_tofd2_diMuons_BR, 2) # total BR for process
    
    return ZD_tofd1fd1bar_BR, fd1_tofd2_diMuons_BR, processBR

def getAllBR(fullPath):
    #####################################################
    # Function that retrieves all BR information for ZD #
    # and fd1 from the parameter card (param_card.dat). #
    # Accepts                                           #
    # no args.                                          #
    # Returns:                                          #
    # no values                                         #
    #####################################################
    fileName = baseMG5Dir + "/" + modelDir + "/param_card.dat"
    with open(fileName, 'r') as file:
        contents = file.readlines()
    
    ZDMass  = int(float(contents[46].split()[1]))
    print("ZD mass: {}".format(ZDMass))
    fd1Mass = int(float(contents[47].split()[1]))
    print("fd1 mass: {}".format(fd1Mass))
    fd2Mass = int(float(contents[48].split()[1]))
    print("fd2 mass: {}".format(fd2Mass))

    resultsFileName = "allBranchingRatios_MZD-%d_mfd1-%d_mfd2-%d.txt" % (ZDMass, fd1Mass, fd2Mass) 
    os.system("cp " + fileName + " " + baseMG5Dir + "/" + fullPath + "/" + resultsFileName )

    '''
    # BRs for ZD decay
    ZD_to_fd1fd2         = float(contents[173].split()[0])
    ZD_to_fd1fd1bar      = float(contents[174].split()[0])
    ZD_to_u_ubar         = float(contents[175].split()[0])
    ZD_to_c_cbar         = float(contents[176].split()[0])
    ZD_to_e_ebar         = float(contents[177].split()[0])
    ZD_to_mu_mubar       = float(contents[178].split()[0])
    ZD_to_tau_taubar     = float(contents[179].split()[0])
    ZD_to_d_dbar         = float(contents[180].split()[0])
    ZD_to_s_sbar         = float(contents[181].split()[0])
    ZD_to_b_bbar         = float(contents[182].split()[0])
    ZD_to_nutau_nutaubar = float(contents[183].split()[0])
    ZD_to_numu_numubar   = float(contents[184].split()[0])
    ZD_to_nue_nuebar     = float(contents[185].split()[0])
    print("ZD to nue nuebar: {}".format(ZD_to_nue_nuebar))


    # BRs for fd1 decay
    fd1_to_u_ubar         = float(contents[190].split()[0])
    fd1_to_e_ebar         = float(contents[191].split()[0])
    fd1_to_mu_mubar       = float(contents[192].split()[0])
    fd1_to_d_dbar         = float(contents[193].split()[0])
    fd1_to_s_sbar         = float(contents[194].split()[0])
    fd1_to_nue_nuebar     = float(contents[195].split()[0])
    fd1_to_numu_numubar   = float(contents[196].split()[0])
    fd1_to_nutau_nutaubar = float(contents[197].split()[0])
    fd1_to_c_cbar         = float(contents[198].split()[0])
    
    # BRs for ZD decay
    ZD_to_fd1fd2         = contents[173].split()[0]
    ZD_to_fd1fd1bar      = contents[174].split()[0]
    ZD_to_u_ubar         = contents[175].split()[0]
    ZD_to_c_cbar         = contents[176].split()[0]
    ZD_to_e_ebar         = contents[177].split()[0]
    ZD_to_mu_mubar       = contents[178].split()[0]
    ZD_to_tau_taubar     = contents[179].split()[0]
    ZD_to_d_dbar         = contents[180].split()[0]
    ZD_to_s_sbar         = contents[181].split()[0]
    ZD_to_b_bbar         = contents[182].split()[0]
    ZD_to_nutau_nutaubar = contents[183].split()[0]
    ZD_to_numu_numubar   = contents[184].split()[0]
    ZD_to_nue_nuebar     = contents[185].split()[0]
    print("ZD to nue nuebar: {}".format(ZD_to_nue_nuebar))


    # BRs for fd1 decay
    fd1_to_u_ubar         = contents[190].split()[0]
    fd1_to_e_ebar         = contents[191].split()[0]
    fd1_to_mu_mubar       = contents[192].split()[0]
    fd1_to_d_dbar         = contents[193].split()[0]
    fd1_to_s_sbar         = contents[194].split()[0]
    fd1_to_nue_nuebar     = contents[195].split()[0]
    fd1_to_numu_numubar   = contents[196].split()[0]
    fd1_to_nutau_nutaubar = contents[197].split()[0]
    fd1_to_c_cbar         = contents[198].split()[0]

    print("fd1 to c cbar: {}".format(fd1_to_c_cbar))

    resultsFileName = "allBranchingRatios_MZD-%d_mfd1-%d_mfd2-%d.txt" % (ZDMass, fd1Mass, fd2Mass) # file for all BRs
    with open(baseMG5Dir + "/" + fullPath + "/" + resultsFileName, 'w') as file:
        file.write("%s    %s    %s    %s    %s    %s    %s    %s    %s    %s    %s    %s    %s    %s    %s    %s    %s    %s    %s    %s    %s    %s\n" 
            % ( "ZD_to_fd1fd2", "ZD_to_fd1fd1bar", "ZD_to_u_ubar", "ZD_to_c_cbar", "ZD_to_e_ebar", "ZD_to_mu_mubar", "ZD_to_tau_taubar", "ZD_to_d_dbar", "ZD_to_s_sbar",
                "ZD_to_b_bbar", "ZD_to_nutau_nutaubar", "ZD_to_numu_numubar", "ZD_to_nue_nuebar", "fd1_to_u_ubar", "fd1_to_e_ebar", "fd1_to_mu_mubar", "fd1_to_d_dbar",
                "fd1_to_s_sbar", "fd1_to_nue_nuebar", "fd1_to_numu_numubar", "fd1_to_nutau_nutaubar", "fd1_to_c_cbar" ))     

        file.write("%s    %s    %s    %s    %s    %s    %s    %s    %s    %s    %s    %s    %s    %s    %s    %s    %s    %s    %s    %s    %s    %s" 
            % ( ZD_to_fd1fd2, ZD_to_fd1fd1bar, ZD_to_u_ubar, ZD_to_c_cbar, ZD_to_e_ebar, ZD_to_mu_mubar, ZD_to_tau_taubar, ZD_to_d_dbar, ZD_to_s_sbar,
                ZD_to_b_bbar, ZD_to_nutau_nutaubar, ZD_to_numu_numubar, ZD_to_nue_nuebar, fd1_to_u_ubar, fd1_to_e_ebar, fd1_to_mu_mubar, fd1_to_d_dbar,
                fd1_to_s_sbar, fd1_to_nue_nuebar, fd1_to_numu_numubar, fd1_to_nutau_nutaubar, fd1_to_c_cbar))
    '''

def getAcc(targetDir):
    #####################################################
    # Function that retrieves the geometric and         #
    # kinematic acceptances. Accepts no arguments.      #
    # from the parameter card (param_card.dat). Accepts #
    # no args.                                          #
    # Returns:                                          #
    # totalEff - Float; cumulative efficiency           #
    #            (all cuts)                             #
    #####################################################
    fileName = baseMG5Dir + "/"+ ma5Dir + "/" + targetDir + "/" + ma5_latex + "/main.tex"

    with open(fileName, 'r') as file:
        contents = file.readlines()

    totalEff = float(contents[263].split()[-3]) # gets cumulative efficiency at cut 8

    return totalEff

def prodXsecCalc(fullPath, mass, verbose, quiet):
    #####################################################
    # Function that performs the prod. xsec             #
    # calculation. Returns nothing.                     #
    # Arguments:                                        #
    # fullPath - String; (relative) path to the results #
    #            dir                                    #
    # mass     - Float; ZD mass                         #
    # verbose  - Boolean; verbose mode                  #            
    #####################################################
    fileXsec = open(fullPath + "/prod_crossX.txt", "a") # open file in append mode
    if verbose == True:
        print(52 * "*")
        print(13 * "*" + " Preparing parameters file " + 13 * "*")
        print(52 * "*")
        print("Editing param_card.dat with ZD mass: %d" % mass)
        editParamCard(mass, default_fd1Mass, default_fd2Mass, modelDir)
    else:
        editParamCard(mass, default_fd1Mass, default_fd2Mass, modelDir)
        
    if verbose == True:
        print(52 * "*")
        print(14 * "*" + " Preparing process card " + 14 * "*")
        print(52 * "*")
        print("Editing proc_card.dat for process %s" % processes[0])
        editProcCard(processes[0], computeWidths = False)
        print("Process card successfully updated")
        
        print(52 * "*")
        print(12 * "*" +  " Generating process %s " % processes[0] + 11 * "*")          
        print(52 * "*")
        
        os.system('./bin/mg5_aMC proc_card.dat')
        
        print(52 * "*")
        print(10 * "*" + " Generating events for %s " % processes[0] + 10 * "*")
        print(52 * "*")
        
        os.system('./' + outputDir + '/bin/generate_events')
        
        print(52 * "*")
        print(4 * "*" + " Getting number of events and cross-section " + 4 * "*")
        print(52 * "*")
        
        numEvents, prodXsec = getXsec(quiet)
        fileXsec.write("%f    %d    %f\n" % (mass, numEvents, prodXsec))
        fileXsec.close()
    
    else:
        editProcCard(processes[0])
        if quiet == True:
            os.system("./bin/mg5_aMC proc_card.dat" + suppress) # generate process
            os.system('./' + outputDir + '/bin/generate_events' + suppress) # generate events
            numEvents, prodXsec = getXsec(quiet)
            fileXsec.write("%f    %d    %f\n" % (mass, numEvents, prodXsec))
            fileXsec.close()
        else:
            os.system('./bin/mg5_aMC proc_card.dat') # generate process
            os.system('./' + outputDir + '/bin/generate_events') # generate events
            numEvents, prodXsec = getXsec(quiet)
            fileXsec.write("%f    %d    %f\n" % (mass, numEvents, prodXsec))
            fileXsec.close()

def accCalc(fullPath, verbose, quiet, mass, fd1Mass, fd2Mass):
    #####################################################
    # Function that performs the acceptances            #
    # calculation. Returns nothing.                     #
    # Arguments:                                        #
    # fullPath - String; (relative) path to the results #
    #            dir                                    #
    # verbose  - Boolean; verbose mode                  #
    # mass     - Float; ZD mass                         #
    # fd1Mass  - Float; fd1 mass                        #
    # fd2Mass  - Float; fd2 mass                        #      
    #####################################################
    if verbose == True:
        fileAcc = open(fullPath + "/acceptances.txt", "a")

        print(52 * "*")
        print(13 * "*" + " Generating events " + 13 * "*")
        print(52 * "*")

        os.system('./' + outputDir + '/bin/generate_events')

        print(52 * "*")
        print(3 * "*" + " Preapring  " + 3 * "*")
        print(52 * "*")

        os.system("cp " + outputEventDir + "/unweighted_events.lhe.gz " + ma5_lheDir)

        os.chdir(ma5Dir) # change to ma5 base dir 

        print(52 * "*")
        print(13 * "*" + " Applying cuts " + 13 * "*")
        print(52 * "*")

        os.system("./bin/ma5 -s " + ma5_script)

        print(52 * "*")
        print(4 * "*" + " Getting geometric and kinematic acceptances " + 4 * "*")
        print(52 * "*")

        ma5_outputDir = getRecentDir()

        totalAcc = getAcc(ma5_outputDir)
        fileAcc.write("%f    %f    %f    %f\n" % (mass, fd1Mass, fd2Mass, totalAcc))
        fileAcc.close()

        os.chdir(baseMG5Dir)

    else:
        if quiet == True:
            #print("opening file: " + fullPath + "/acceptances.txt")
            fileAcc = open(fullPath + "/acceptances.txt", "a")
            #print("generating events")
            os.system('./' + outputDir + '/bin/generate_events' + suppress)
            #print("copying lhe file to ma5 directory")
            os.system("cp " + outputEventDir + "/unweighted_events.lhe.gz " + ma5_lheDir + suppress)
            #print("Changing dirs to ma5 dir")
            os.chdir(ma5Dir) # change to ma5 base dir 
            #print("running ma5 analysis")
            os.system("./bin/ma5 -s " + ma5_script + suppress)
            
            ma5_outputDir = getRecentDir()
            print("Recent dir: " + ma5_outputDir)
            totalAcc = getAcc(ma5_outputDir)
            print("total acceptance: %f" % totalAcc )
            fileAcc.write("%f    %f    %f    %f\n" % (mass, fd1Mass, fd2Mass, totalAcc))
            fileAcc.close()

            os.chdir(baseMG5Dir)
        else:
            fileAcc = open(fullPath + "/acceptances.txt", "a")
            os.system('./' + outputDir + '/bin/generate_events' + suppress)
            os.system("cp " + outputEventDir + "/unweighted_events.lhe.gz " + ma5_lheDir + suppress)
            os.chdir(ma5Dir) # change to ma5 base dir 
            os.system("./bin/ma5 -s " + ma5_script + suppress)
            ma5_outputDir = getRecentDir()
            totalAcc = getAcc(ma5_outputDir)
            fileAcc.write("%f    %f    %f    %f\n" % (mass, fd1Mass, fd2Mass, totalAcc))
            fileAcc.close()


def branchingRatioCalc(fullPath, mass, processes, acceptances, allBR, fd1MassScan, fd2MassScan, verbose, quiet, massList_fd1 = None, massList_fd2 = None):
    #####################################################
    # Function that performs the BR calculation.        #
    # Returns nothing.                                  #
    # Arguments:                                        #
    # mass      - Float; ZD mass                        #
    # processes   - List of strings; List of processes  #
    # acceptances - Boolean; perform acceptance         #
    #               calculation                         #
    # allBR       - Boolean; retrieve all BRs for ZD    #
    #             - and fd1                             #
    # fd1MassScan - Boolean; perform fd2 mass scan      #
    # fd2MassScan - Boolean; perform fd2 mass scan      #
    # verbose     - Boolean; verbose mode               #
    # Optional arguments:                               #
    # massList_fd1 - Numpy array (1D); Array of fd1     #
    #                masses                             #
    # massList_fd2 - Numpy array (1D); Array of fd1     #
    #                masses                             #    
    #####################################################
    if fd1MassScan == False and fd2MassScan == False:
        fileBR = open(fullPath + "/branching_ratio.txt", "a")

        if verbose == True:
            print(52 * "*")
            print(14 * "*" + " Preparing process card " + 14 * "*")
            print(52 * "*")
            print("Editing proc_card.dat for process %s" % processes[1])
            editProcCard(processes[1])
            print("Process card successfully updated")
            
            print(72 * "*")
            print("* Generating process %s *" % processes[1])            
            print(72 * "*")
        
            os.system('./bin/mg5_aMC proc_card.dat')
            
            print(74 * "*")
            print(26 * "*" + "Determining branching ratio" + " " + 26 * "*")           
            print(26 * "*" + " for %s " % processes[1] + " " + 26 * "*")
            
            ZD_tofd1fd1bar_BR, fd1_tofd2_diMuons_BR, processBR = getBR()
            
            fileBR.write("%f    %f    %f    %f    %f    %f\n" % (mass, default_fd1Mass, default_fd2Mass, ZD_tofd1fd1bar_BR, fd1_tofd2_diMuons_BR, processBR))
            fileBR.close()

            if allBR == True:
                print(74 * "*")
                print(26 * "*" + "Retrieving all branching ratios" + " " + 26 * "*")           
                print(26 * "*" + " for %s " % processes[1] + " " + 26 * "*")
                getAllBR(fullPath)
            else:
                pass

            if acceptances == True:
                accCalc(fullPath, verbose, mass, default_fd1Mass, default_fd2Mass)
            else:
                pass

        else:
            if quiet == True:
                editProcCard(processes[1], computeWidths = True)
                os.system('./bin/mg5_aMC proc_card.dat' + suppress)
                ZD_tofd1fd1bar_BR, fd1_tofd2_diMuons_BR, processBR = getBR()
                fileBR.write("%f    %f    %f    %f    %f    %f\n" % (mass, default_fd1Mass, default_fd2Mass, ZD_tofd1fd1bar_BR, fd1_tofd2_diMuons_BR, processBR))
                fileBR.close()

                if allBR == True:
                    getAllBR(fullPath)
                else:
                    pass

                if acceptances == True:
                    accCalc(fullPath, verbose, mass, quiet, default_fd1Mass, default_fd2Mass)
                else:
                    pass
            else:
                editProcCard(processes[1], computeWidths = True)
                os.system('./bin/mg5_aMC proc_card.dat')
                ZD_tofd1fd1bar_BR, fd1_tofd2_diMuons_BR, processBR = getBR()
                fileBR.write("%f    %f    %f    %f    %f    %f\n" % (mass, default_fd1Mass, default_fd2Mass, ZD_tofd1fd1bar_BR, fd1_tofd2_diMuons_BR, processBR))
                fileBR.close()

                if acceptances == True:
                    accCalc(fullPath, verbose, mass, quiet, default_fd1Mass, default_fd2Mass)
                else:
                    pass

                if allBR == True:
                    getAllBR(fullPath)
                else:
                    pass

    elif fd1MassScan == True and fd2MassScan == False:
        counter_fd1 = 0
        if verbose == True:
            for fd1Mass in massList_fd1:
                if 2 * fd1Mass > mass: # skip kinematically impossible fd1 masses
                    continue

                fileBR = open(fullPath + "/branching_ratio.txt", "a")
                print(52 * "*")
                print(13 * "*" + " Preparing parameters file " + 13 * "*")
                print(52 * "*")
                print("Editing param_card.dat with fd1 mass: %d" % fd1Mass)
                editParamCard(mass, fd1Mass, default_fd2Mass, modelDir)
                
                print(52 * "*")
                print(14 * "*" + " Preparing process card " + 14 * "*")
                print(52 * "*")
                print("Editing proc_card.dat for process %s" % processes[1])
                editProcCard(processes[1], computeWidths = True)
                print("Process card successfully updated")
                
                print(72 * "*")
                print("* Generating process %s *" % processes[1])            
                print(72 * "*")
            
                os.system('./bin/mg5_aMC proc_card.dat')
                
                print(74 * "*")
                print(26 * "*" + "Determining branching ratio" + " " + 26 * "*")           
                print(26 * "*" + " for %s " % processes[1] + " " + 26 * "*")
                
                ZD_tofd1fd1bar_BR, fd1_tofd2_diMuons_BR, processBR = getBR()
                
                fileBR.write("%f    %f    %f    %f    %f    %f\n" % (mass, fd1Mass, default_fd2Mass, ZD_tofd1fd1bar_BR, fd1_tofd2_diMuons_BR, processBR))
                fileBR.close()

                if allBR == True:
                    print(74 * "*")
                    print(26 * "*" + "Retrieving all branching ratios" + " " + 26 * "*")           
                    print(26 * "*" + " for %s " % processes[1] + " " + 26 * "*")
                    getAllBR(fullPath)
                else:
                    pass

                if acceptances == True:
                    accCalc(fullPath, verbose, quiet, mass, fd1Mass, default_fd2Mass)
                else:
                    pass

        else:
            if quiet == True:
                for fd1Mass in tqdm(massList_fd1, ascii = True, desc = "Scanning over fd1 mass"):
                    print('\n')
                    if 2 * fd1Mass > mass: # skip kinematically impossible fd1 masses
                        continue

                    fileBR = open(fullPath + "/branching_ratio.txt", "a")
                    editParamCard(mass, fd1Mass, default_fd2Mass, modelDir)
                    editProcCard(processes[1], computeWidths = True)
                    os.system('./bin/mg5_aMC proc_card.dat' + suppress)
                    ZD_tofd1fd1bar_BR, fd1_tofd2_diMuons_BR, processBR = getBR()
                    fileBR.write("%f    %f    %f    %f    %f    %f\n" % (mass, fd1Mass, default_fd2Mass, ZD_tofd1fd1bar_BR, fd1_tofd2_diMuons_BR, processBR))
                    fileBR.close()

                    if acceptances == True:
                        accCalc(fullPath, verbose, quiet, mass, default_fd1Mass, default_fd2Mass)
                    else:
                        pass
            else:
                for fd1Mass in massList_fd1:
                    if 2 * fd1Mass > mass: # skip kinematically impossible fd1 masses
                        continue

                    fileBR = open(fullPath + "/branching_ratio.txt", "a")
                    editParamCard(mass, fd1Mass, default_fd2Mass, modelDir)
                    editProcCard(processes[1], computeWidths = True)
                    os.system('./bin/mg5_aMC proc_card.dat')
                    ZD_tofd1fd1bar_BR, fd1_tofd2_diMuons_BR, processBR = getBR()
                    fileBR.write("%f    %f    %f    %f    %f    %f\n" % (mass, fd1Mass, default_fd2Mass, ZD_tofd1fd1bar_BR, fd1_tofd2_diMuons_BR, processBR))
                    fileBR.close()

                    if acceptances == True:
                        accCalc(fullPath, verbose, quiet, mass, default_fd1Mass, default_fd2Mass)
                    else:
                        pass

    elif fd1MassScan == False and fd2MassScan == True:
        if verbose == True:
            for fd2Mass in massList_fd2:
                if fd2Mass + (2 * muonMass) > default_fd1Mass: # skip kinematically impossible fd2 masses
                    continue

                fileBR = open(fullPath + "/branching_ratio.txt", "a")

                print(52 * "*")
                print(13 * "*" + " Preparing parameters file " + 13 * "*")
                print(52 * "*")
                print("Editing param_card.dat with fd2 mass: %d" % fd1Mass)
                editParamCard(mass, default_fd1Mass, fd2Mass, modelDir)

                print(52 * "*")
                print(14 * "*" + " Preparing process card " + 14 * "*")
                print(52 * "*")
                print("Editing proc_card.dat for process %s" % processes[1])
                editProcCard(processes[1], computeWidths = True)
                print("Process card successfully updated")
                
                print(72 * "*")
                print("* Generating process %s *" % processes[1])            
                print(72 * "*")
            
                os.system('./bin/mg5_aMC proc_card.dat')
                
                print(74 * "*")
                print(26 * "*" + "Determining branching ratio" + " " + 26 * "*")           
                print(26 * "*" + " for %s " % processes[1] + " " + 26 * "*")
                
                ZD_tofd1fd1bar_BR, fd1_tofd2_diMuons_BR, processBR = getBR()
                
                fileBR.write("%f    %f    %f    %f    %f    %f\n" % (mass, default_fd1Mass, fd2Mass, ZD_tofd1fd1bar_BR, fd1_tofd2_diMuons_BR, processBR))
                fileBR.close()

                if allBR == True:
                    print(74 * "*")
                    print(26 * "*" + "Retrieving all branching ratios" + " " + 26 * "*")           
                    print(26 * "*" + " for %s " % processes[1] + " " + 26 * "*")
                    getAllBR(fullPath)
                else:
                    pass

                if acceptances == True:
                    accCalc(fullPath, verbose, quiet, mass, default_fd1Mass, fd2Mass)
                else:
                    pass
        else:
            if quiet == True:
                for fd2Mass in tqdm(massList_fd2, ascii = True, desc = "Scanning over fd2 mass"):
                    print('\n')
                    if fd2Mass + (2 * muonMass) > default_fd1Mass: # skip kinematically impossible fd2 masses
                        continue

                    fileBR = open(fullPath + "/branching_ratio.txt", "a")
                    editParamCard(mass, default_fd1Mass, fd2Mass, modelDir)
                    editProcCard(processes[1], computeWidths = True)
                    os.system('./bin/mg5_aMC proc_card.dat' + suppress)
                    ZD_tofd1fd1bar_BR, fd1_tofd2_diMuons_BR, processBR = getBR()
                    fileBR.write("%f    %f    %f    %f    %f    %f\n" % (mass, fd1Mass, default_fd2Mass, ZD_tofd1fd1bar_BR, fd1_tofd2_diMuons_BR, processBR))
                    fileBR.close()

                    if allBR == True:
                        getAllBR(fullPath)
                    else:
                        pass

                    if acceptances == True:
                        accCalc(fullPath, verbose, quiet, mass, default_fd1Mass, fd2Mass)
                    else:
                        pass
            else:
                for fd2Mass in massList_fd2:
                    if fd2Mass + (2 * muonMass) > default_fd1Mass: # skip kinematically impossible fd2 masses
                        continue

                    fileBR = open(fullPath + "/branching_ratio.txt", "a")
                    editParamCard(mass, default_fd1Mass, fd2Mass, modelDir)
                    editProcCard(processes[1], computeWidths = True)
                    os.system('./bin/mg5_aMC proc_card.dat')
                    ZD_tofd1fd1bar_BR, fd1_tofd2_diMuons_BR, processBR = getBR()
                    fileBR.write("%f    %f    %f    %f    %f    %f\n" % (mass, fd1Mass, default_fd2Mass, ZD_tofd1fd1bar_BR, fd1_tofd2_diMuons_BR, processBR))
                    fileBR.close()

                    if allBR == True:
                        getAllBR(fullPath)
                    else:
                        pass

                    if acceptances == True:
                        accCalc(fullPath, verbose, quiet, mass, default_fd1Mass, fd2Mass)
                    else:
                        pass   

    elif fd1MassScan == True and fd2MassScan == True:
        if verbose == True:
            for fd1Mass in massList_fd1:
                if 2 * fd1Mass > mass: # skip kinematically impossible fd1 masses
                    continue

                for fd2Mass in massList_fd2:

                    if fd1Mass < fd2Mass + (2 * muonMass): # skip the rest of the loop if fd1 < fd2 + dimuon mass (kinematically impossible)
                        continue

                    fileBR = open(fullPath + "/branching_ratio.txt", "a")

                    print(52 * "*")
                    print(13 * "*" + " Preparing parameters file " + 13 * "*")
                    print(52 * "*")
                    print("Editing param_card.dat with fd2 mass: %d" % fd1Mass)
                    editParamCard(mass, fd1Mass, fd2Mass, modelDir)

                    print(52 * "*")
                    print(14 * "*" + " Preparing process card " + 14 * "*")
                    print(52 * "*")
                    print("Editing proc_card.dat for process %s" % processes[1])
                    editProcCard(processes[1], computeWidths = True)
                    print("Process card successfully updated")
                    
                    print(72 * "*")
                    print("* Generating process %s *" % processes[1])            
                    print(72 * "*")
                
                    os.system('./bin/mg5_aMC proc_card.dat')
                    
                    print(74 * "*")
                    print(26 * "*" + "Determining branching ratio" + " " + 26 * "*")           
                    print(26 * "*" + " for %s " % processes[1] + " " + 26 * "*")
                    
                    ZD_tofd1fd1bar_BR, fd1_tofd2_diMuons_BR, processBR = getBR()
                    
                    fileBR.write("%f    %f    %f    %f    %f    %f\n" % (mass, fd1Mass, fd2Mass, ZD_tofd1fd1bar_BR, fd1_tofd2_diMuons_BR, processBR))
                    fileBR.close()

                    if allBR == True:
                        print(74 * "*")
                        print(26 * "*" + "Retrieving all branching ratios" + " " + 26 * "*")           
                        print(26 * "*" + " for %s " % processes[1] + " " + 26 * "*")
                        getAllBR(fullPath)
                    else:
                        pass

                    if acceptances == True:
                        accCalc(fullPath, verbose, quiet, mass, fd1Mass, fd2Mass)
                    else:
                        pass
        else:
            if quiet == True:
                for fd1Mass in tqdm(massList_fd1, ascii = True, desc = "Scanning over fd1 masses"):
                    print('\n')
                    if 2 * fd1Mass > mass: # skip kinematically impossible fd1 masses
                        continue

                    for fd2Mass in tqdm(massList_fd2, ascii = True, desc = "Scanning over fd2 masses"):
                        print('\n')
                        if fd1Mass < fd2Mass + (2 * muonMass): # skip the rest of the loop if fd1 < fd2 + dimuon mass (kinematically impossible)
                            continue

                        fileBR = open(fullPath + "/branching_ratio.txt", "a")
                        editParamCard(mass, fd1Mass, fd2Mass, modelDir)
                        editProcCard(processes[1], computeWidths = True)
                        os.system('./bin/mg5_aMC proc_card.dat' + suppress)
                        ZD_tofd1fd1bar_BR, fd1_tofd2_diMuons_BR, processBR = getBR()           
                        fileBR.write("%f    %f    %f    %f    %f    %f\n" % (mass, fd1Mass, fd2Mass, ZD_tofd1fd1bar_BR, fd1_tofd2_diMuons_BR, processBR))
                        fileBR.close()

                        if allBR == True:
                            getAllBR(fullPath)
                        else:
                            pass

                        if acceptances == True:
                            accCalc(fullPath, verbose, quiet, mass, fd1Mass, fd2Mass)
                        else:
                            pass
            else:
                for fd1Mass in massList_fd1:
                    if 2 * fd1Mass > mass: # skip kinematically impossible fd1 masses
                        continue

                    for fd2Mass in massList_fd2:
                        if fd1Mass < fd2Mass + (2 * muonMass): # skip the rest of the loop if fd1 < fd2 + dimuon mass (kinematically impossible)
                            continue

                        fileBR = open(fullPath + "/branching_ratio.txt", "a")
                        editParamCard(mass, fd1Mass, fd2Mass, modelDir)
                        editProcCard(processes[1], computeWidths = True)
                        os.system('./bin/mg5_aMC proc_card.dat')
                        ZD_tofd1fd1bar_BR, fd1_tofd2_diMuons_BR, processBR = getBR()           
                        fileBR.write("%f    %f    %f    %f    %f    %f\n" % (mass, fd1Mass, fd2Mass, ZD_tofd1fd1bar_BR, fd1_tofd2_diMuons_BR, processBR))
                        fileBR.close()

                        if allBR == True:
                            getAllBR(fullPath)
                        else:
                            pass

                        if acceptances == True:
                            accCalc(fullPath, verbose, quiet, mass, fd1Mass, fd2Mass)
                        else:
                            pass    

def runScans(fullPath, massList, prodCrossSection, branchingRatio, acceptances, allBR, verbose, quiet, fd1MassScan, fd2MassScan, massList_fd1 = None, massList_fd2 = None):
    #####################################################
    # Function that performs the BR calculation.        #
    # Returns nothing.                                  #
    # Arguments:                                        #
    # fullPath         - String; relative path to the   #
    #                    results directory              #
    # massList         - Numpy array (1D); array of ZD  #
    #                    masses to scan                 #    
    # prodCrossSection - Boolean; perform prod. xsec    #
    #                    calculation                    #
    # branchingRatio   - Boolean; perform BR            #
    #                    calculation                    #
    # acceptances      - Boolean; perform geo. + kin.   #
    #                    acceptance calculation         #
    # allBR            - Retrieve all BRs for ZD and    #
    #                    fd1                            #
    # verbose          - Boolean; verbose mode          #
    # fd1MassScan      - Boolean; perform fd2 mass scan #
    # fd2MassScan      - Boolean; perform fd2 mass scan #
    # Optional arguments:                               #
    # massList_fd1 - Numpy array (1D); Array of fd1     #
    #                masses.                            #
    # massList_fd2 - Numpy array (1D); Array of fd1     #
    #                masses                             #    
    #####################################################    
    prepareDataFiles(fullPath, prodCrossSection, branchingRatio, acceptances) # prepare results file for cross section and BR

    for mass in tqdm(massList, ascii = True, desc = "Scanning over ZD mass"):

        if prodCrossSection == True:
            print('\n')
            prodXsecCalc(fullPath, mass, verbose, quiet)
        else:
            pass
            
        if branchingRatio == True:
            print('\n')
            branchingRatioCalc(fullPath, mass, processes, acceptances, allBR, fd1MassScan, fd2MassScan, verbose, quiet, massList_fd1, massList_fd2)
        else:
            pass
    
if __name__ == "__main__":

    # arg parser
    parser = argparse.ArgumentParser(description="Script to get BR, production cross-section, and partial cross sections for various MZD for fD model")
    parser.add_argument("-ml", "--massLow", action="store", dest="massLow", help="Lower bound of MZD to scan (only integer value)")
    parser.add_argument("-mh", "--massHigh", action="store", dest="massHigh", help="Upper bound of MZD to scan (only integer value)")
    parser.add_argument("-mi", "--massIncrement", action="store", dest="massInc", help="MZD increment (only integer value)")
    parser.add_argument("-n", "--name", action="store", dest="runName", help="Name of results directory")
    parser.add_argument("-q", "--quiet", action="store_true", dest="quiet", default = False, help="Quiet mode: suppress output from MG5 and MA5. Also prints progress bar.")    
    parser.add_argument("-v", "--verbose", action="store_true", dest="verbose", default = False, help="verbose mode")
    parser.add_argument("-b", "--branchingRatio", action="store_true", dest="branchingRatio", default = False, help="Calculate branching ratios")
    parser.add_argument("-p", "--prodCrossSection", action="store_true", dest="prodCrossSection", default = False, help="Calculate production cross-section")
    parser.add_argument("-a", "--acceptances", action="store_true", dest="acceptances", default = False, help="Calculate production cross-section")
    parser.add_argument("-f1l", "--fd1massLow", action="store", dest="fd1massLow", help="Lower bound of fd1 mass to scan (only integer value)")
    parser.add_argument("-ab", "--allBranchingRatios", action="store_true", dest="allBR", default = False, help="Retrieve all branching ratios")
    parser.add_argument("-f1h", "--fd1massHigh", action="store", dest="fd1massHigh", help="Upper bound of fd1 mass to scan (only integer value)")
    parser.add_argument("-f1i", "--fd1massInc", action="store", dest="fd1massInc", help="fd1 mass increment (only integer value)")
    parser.add_argument("-f2l", "--fd2massLow", action="store", dest="fd2massLow", help="Lower bound of fd2 mass to scan (only integer value)")
    parser.add_argument("-f2h", "--fd2massHigh", action="store", dest="fd2massHigh", help="Upper bound of fd2 mass to scan (only integer value)")
    parser.add_argument("-f2i", "--fd2massInc", action="store", dest="fd2massInc", help="fd2 mass increment (only integer value)")

    args = parser.parse_args()    
    
    # Check arguments
    if os.path.isdir(massScanDir + "/" + args.runName):
        proceed = str(input(colors.RED + "Directory " + args.runName + " already exists. Do you wish to overwrite the contents? (y/n): " + colors.ENDC))
        if proceed == "y":
            doubleCheck = str(input(colors.YELLOW + "You selected yes. Are you sure? This will erase any existing data files. (y/n): " + colors.ENDC))
            if doubleCheck == "y":
                print(colors.GREEN + "Directory " + args.runName + " initialized. Starting calculations now." + colors.ENDC)
            else:
                print(colors.YELLOW + "Exiting" + colors.ENDC)
                sys.exit()
        else:
            print(colors.YELLOW + "Exiting" + colors.ENDC)
            sys.exit()

    if args.massLow is None or args.massHigh is None:
        print(colors.YELLOW + "Must specify the mass range to scan" + colors.ENDC)
        print(colors.YELLOW + "Exiting" + colors.ENDC)
        sys.exit()

    if args.massInc is None:
        print(colors.YELLOW + "Must specify the mass increment" + colors.ENDC)
        print(colors.YELLOW + "Exiting" + colors.ENDC)
        sys.exit()

    if args.fd1massLow is not None and args.fd1massLow is None and args.fd1massInc is None:
        print (colors.YELLOW + "If scanning over fd1 mass the mass range must be specified" + colors.ENDC)
        print(colors.YELLOW + "Exiting" + colors.ENDC)
        sys.exit()

    if args.fd2massLow is not None and args.fd2massLow is None and args.fd2massInc is None:
        print (colors.YELLOW + "If scanning over fd2 mass the mass range must be specified" + colors.ENDC)
        print(colors.YELLOW + "Exiting" + colors.ENDC)
        sys.exit()

    if args.branchingRatio and args.prodCrossSection is None:
        print (colors.YELLOW + "Must select branching ratio and/or production cross-section" + colors.ENDC)
        print(colors.YELLOW + "Exiting" + colors.ENDC)
        sys.exit()

    if args.runName is None:
        print (colors.YELLOW + "Must specify run name" + colors.ENDC)
        sys.exit()

    if args.acceptances == True and args.branchingRatio == False:
        print (colors.YELLOW + "Must perform the branching ratio calculation to calculate kinematic and geometric acceptances" + colors.ENDC)
        print(colors.YELLOW + "Exiting" + colors.ENDC)
        sys.exit()

    if args.allBR == True and args.branchingRatio == False:
        print (colors.YELLOW + "Must perform the branching ratio calculation to retrieve all branching ratios" + colors.ENDC)
        print(colors.YELLOW + "Exiting" + colors.ENDC)
        sys.exit()

    if args.quiet == True and args.verbose == True:
        print (colors.YELLOW + "Cannot select quiet and verbose mode at the same time" + colors.ENDC)
        print(colors.YELLOW + "Exiting" + colors.ENDC)
        sys.exit()

    if args.fd1massLow is not None and args.fd1massHigh is not None and args.fd1massInc is not None:
        fd1MassScan  = True
        massList_fd1 = np.arange(float(args.fd1massLow), float(args.fd1massHigh), float(args.fd1massInc), dtype = float) # initialize fd1 mass list
    else:
        fd1MassScan  = False
        massList_fd1 = None
    if args.fd2massLow is not None and args.fd2massHigh is not None and args.fd2massInc is not None:
        fd2MassScan = True 
        massList_fd2 = np.arange(float(args.fd2massLow), float(args.fd2massHigh), float(args.fd2massInc), dtype = float) # initialize fd1 mass list
    else:
        fd2MassScan  = False
        massList_fd2 = None
        
    massList = np.arange(float(args.massLow), float(args.massHigh), float(args.massInc), dtype = float) # initialize ZD mass list

    full_path = newScanSetup(args.runName)

    try:
        runScans(full_path, massList, args.prodCrossSection, args.branchingRatio, args.acceptances, args.allBR, args.verbose, args.quiet, fd1MassScan, fd2MassScan, massList_fd1, massList_fd2)
    except KeyboardInterrupt:
        print(colors.RED + "\nKeyboard interrupt: exiting..." + colors.ENDC)

    

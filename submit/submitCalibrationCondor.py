#!/usr/bin/env python3

import subprocess, time, sys, os
from methods import *
from datetime import datetime

from optparse import OptionParser                                                                                                                   

parser = OptionParser(usage="%prog [options]")    
parser.add_option("-c", "--create", dest="create", action="store_true", default=False, help="Do not submit the jobs, only create the subfolders")
parser.add_option("-l", "--daemon-local", dest="daemonLocal", action="store_true", default=False, help="Do not submit a job to manage the daemon, do it locally")
parser.add_option("--recover-fill", dest="recoverFill", action="store_true", default=False, help="Before moving to the hadd part of the calibration, first try to recover failed fills")
parser.add_option("-t", "--token-file", dest="tokenFile", type="string", default="", help="File needed to renew token (when daemon running locally)")
parser.add_option("--min-efficiency-recover-fill", dest="minEfficiencyToRecoverFill", type="float", default=0.97, help="Tolerance of EcalNtp loss. Require fraction of good EcalNtp above this number to skip recover")
(options, args) = parser.parse_args()
pwd = os.getcwd()

if ContainmentCorrection == '2017reg':
    os.system("cp /eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero2017/zhicaiz/GBRForest_2017/* ../FillEpsilonPlot/data/")

#-------- create folders --------#

workdir = f"{pwd}/{dirname}"
condordir = f"{workdir}/condor_files"
cfgFillPath = f"{workdir}/cfgFile/Fill"
cfgFitPath  = f"{workdir}/cfgFile/Fit"
cfgHaddPath = f"{workdir}/src/hadd"
srcPath = f"{workdir}/src"

print(f"[calib] Creating local folders ({dirname})")
subprocess.run(['mkdir', '-p', workdir])
subprocess.run(['mkdir', '-p', condordir])
for it in range(nIterations):
    subprocess.run(['mkdir', '-p', f"{condordir}/iter_{it}"])
subprocess.run(['mkdir', '-p', f"{workdir}/cfgFile"])
subprocess.run(['mkdir', '-p', f"{workdir}/CRAB_files"])
for it in range(nIterations):
    subprocess.run(['mkdir', '-p', f"{cfgFillPath}/iter_{it}"])
subprocess.run(['mkdir', '-p', cfgFitPath])
subprocess.run(['mkdir', '-p', srcPath])
for it in range(nIterations):
    subprocess.run(['mkdir', '-p', f"{srcPath}/Fill/iter_{it}"])
subprocess.run(['mkdir', '-p', f"{srcPath}/Fit"])
subprocess.run(['mkdir', '-p', f"{srcPath}/hadd"])

print("[calib] Storing parameters.py for future reference")
subprocess.run(['cp', 'parameters.py', workdir])

# Create folders on EOS / PNFS
if isOtherT2 and storageSite == "T2_BE_IIHE" and isCRAB:
    print("[calib] Creating folders on PNFS")
    subprocess.run(['srmmkdir', f"srm://maite.iihe.ac.be:8443{eosPath}/{dirname}"])
else:
    print("[calib] Creating folders on EOS")
    subprocess.run(['mkdir', '-p', f"{eosPath}/{dirname}"])

for it in range(nIterations):
    if isOtherT2 and storageSite == "T2_BE_IIHE" and isCRAB:
        print(f"[calib]  ---  srmmkdir {eosPath}/{dirname}/iter_{it}")
        subprocess.run(['srmmkdir', f"srm://maite.iihe.ac.be:8443{eosPath}/{dirname}/iter_{it}"])
    else:
        print(f"[calib]  ---  mkdir {eosPath}/{dirname}/iter_{it}")
        subprocess.run(['mkdir', '-p', f"{eosPath}/{dirname}/iter_{it}"])

#-------- fill cfg files --------#
if isCRAB:
    print("--------------This inter-calibration will use CRAB: Good Luck!------------------------")

# open list of input files
with open(inputlist_n) as inputlist_f:
    inputlistbase_v = [x for x in inputlist_f if not x.lstrip().startswith('#')]

print(f"[calib] Total number of files to be processed: {len(inputlistbase_v)}")
print("[calib] Creating cfg Files")

ijob = 0
for it in range(nIterations):
    print(f"[calib]  '-- Fill::Iteration {it}'")
    inputlist_v = inputlistbase_v[:]
    ijob = 0

    NrelJob = float(len(inputlist_v)) / float(ijobmax)
    if (int(NrelJob) - NrelJob) < 0.:
        NrelJob = int(NrelJob) + 1
    Nlist = int(NrelJob / nHadd + 1.0)

    haddSrc_n_s = []
    haddSrc_f_s = []

    print(f"[calib]  '-- Hadd::Number of hadd tasks: {Nlist}  ({nHadd} files per task)'")

    haddSrc_final_n_s = f"{srcPath}/hadd/hadd_iter_{it}_final.list"
    haddSrc_final_f_s = open(haddSrc_final_n_s, 'w')
    for num_list in range(Nlist):
        haddSrc_n_s.append(f"{srcPath}/hadd/hadd_iter_{it}_step_{num_list}.list")
        haddSrc_f_s.append(open(haddSrc_n_s[num_list], 'w'))
        fileToAdd_final_n_s = f"{eosPath}/{dirname}/iter_{it}/{NameTag}epsilonPlots_{num_list}.root\n"
        for nj in range(nHadd):
            nEff = num_list*nHadd + nj
            fileToAdd_n_s = f"{eosPath}/{dirname}/iter_{it}/{NameTag}{outputFile}_{nEff}.root\n"
            if nEff < NrelJob:
                haddSrc_f_s[num_list].write(fileToAdd_n_s)
        haddSrc_final_f_s.write(fileToAdd_final_n_s)
        haddSrc_f_s[num_list].close()
    haddSrc_final_f_s.close()

    # create Hadd cfg file
    dest = eosPath + '/' + dirname + '/iter_' + str(it) + '/'
    for num_list in range(Nlist):
        hadd_cfg_n = cfgHaddPath + "/HaddCfg_iter_" + str(it) + "_job_" + str(num_list) + ".sh"
        hadd_cfg_f = open( hadd_cfg_n, 'w')
        HaddOutput = NameTag + "epsilonPlots_" + str(num_list) + ".root"
        printParallelHaddFAST(hadd_cfg_f, HaddOutput, haddSrc_n_s[num_list], dest, pwd, num_list)
        hadd_cfg_f.close()
        subprocess.run(["chmod", "777", hadd_cfg_n], check=True)
    # print Final hadd
    Fhadd_cfg_n = cfgHaddPath + "/Final_HaddCfg_iter_" + str(it) + ".sh"
    Fhadd_cfg_f = open( Fhadd_cfg_n, 'w')
    printFinalHaddRegroup(Fhadd_cfg_f, haddSrc_final_n_s, dest, pwd)
    Fhadd_cfg_f.close()
    # loop over the whole list
    while (len(inputlist_v) > 0):

        # create cfg file
        fill_cfg_n = cfgFillPath + "/iter_" + str(it) + "/fillEps_iter_" + str(it) + "_job_" + str(ijob) + ".py"
        print(f"\tWriting {fill_cfg_n} ...")
        fill_cfg_f = open(fill_cfg_n, 'w')

        # print first part of the cfg file
        printFillCfg1( fill_cfg_f )
        # loop over the names of the input files to be put in a single cfg
        lastline = min(ijobmax,len(inputlist_v)) - 1
        for line in range(min(ijobmax,len(inputlist_v))):
            ntpfile = inputlist_v.pop(0)
            ntpfile = ntpfile.rstrip()
            if ntpfile != '':
                prefixSourceFileToUse = ""
                if prefixSourceFile not in ntpfile:
                    prefixSourceFileToUse = prefixSourceFile
                if(line != lastline):
                    fill_cfg_f.write("        '" + prefixSourceFileToUse + ntpfile + "',\n")
                else:
                    fill_cfg_f.write("        '" + prefixSourceFileToUse + ntpfile + "'\n")

        # print the last part of the cfg file
        if( isCRAB ):
            printFillCfg2( fill_cfg_f, pwd, it , "", ijob )
        else: 
            printFillCfg2( fill_cfg_f, pwd, it , "/tmp/", ijob )
        fill_cfg_f.close()

        # print source file for batch submission of FillEpsilonPlot task
        fillSrc_n = srcPath + "/Fill/iter_" + str(it) + "/submit_iter_" + str(it) + "_job_" + str(ijob) + ".sh"
        fillSrc_f = open( fillSrc_n, 'w')
        source_s = NameTag +outputFile + "_" + str(ijob) + ".root"
        destination_s = eosPath + '/' + dirname + '/iter_' + str(it) + "/" + source_s
        logpathFill = pwd + "/" + dirname + "/log/iter_" + str(it) + "/fillEps_iter_" + str(it) + "_job_" + str(ijob) + ".log"
        printSubmitSrc(fillSrc_f, fill_cfg_n, "/tmp/" + source_s, destination_s , pwd, logpathFill)
        fillSrc_f.close()

        # make the source file executable
        changePermission = subprocess.Popen(['chmod 777 ' + fillSrc_n], stdout=subprocess.PIPE, shell=True);
        debugout = changePermission.communicate()

        ijob = ijob+1

njobs = ijob

#-------- fit cfg files --------#
nEBindependentXtals = 1699 if foldInSuperModule else 61199
nEB = nEBindependentXtals // nFit
if nEBindependentXtals % nFit != 0:
    nEB += 1
nEE = 14647 // nFit
if 14647 % nFit != 0:
    nEE += 1

if Barrel_or_Endcap == "ONLY_ENDCAP": nEB = 0
if Barrel_or_Endcap == "ONLY_BARREL": nEE = 0

print(f'[calib] Splitting Fit Task: {nEB} jobs on EB, {nEE} jobs on EE')
print('[calib] Creating Fit cfg files')

inListB = [nFit*tmp for tmp in range(nEB)]
finListB = [nFit*tmp+(nFit-1) for tmp in range(nEB)]
inListE = [nFit*tmp for tmp in range(nEE)]
finListE = [nFit*tmp+(nFit-1) for tmp in range(nEE)]

for it in range(nIterations):
    print(f"[calib]  '-- Fit::Iteration {it}'")

    if foldInSuperModule:
        fit_cfg_n = cfgFitPath + f"/fitEpsilonPlot_justFoldSM_iter_{it}.py"
        with open(fit_cfg_n, 'w') as fit_cfg_f:
            printFitCfg(fit_cfg_f, it, "/tmp", 0, 0, "Barrel", 0, justDoHistogramFolding=True)

        fitSrc_n = srcPath + f"/Fit/submit_justFoldSM_iter_{it}.sh"
        with open(fitSrc_n, 'w') as fitSrc_f:
            destination_s = f"{eosPath}/{dirname}/iter_{it}/{NameTag}Barrel_{nFit}_{calibMapName}"
            logpath = f"{pwd}/{dirname}/log/fitEpsilonPlot_justFoldSM_iter_{it}.log"
            printSubmitFitSrc(
                fitSrc_f,
                fit_cfg_n,
                f"/tmp/{NameTag}justFoldSM_{calibMapName}",
                destination_s,
                pwd,
                logpath,
                justDoHistogramFolding=True
            )

        subprocess.run(['chmod', '777', fitSrc_n], check=True)

    for nFit in range(nEB):
        fit_cfg_n = cfgFitPath + f"/fitEpsilonPlot_EB_{nFit}_iter_{it}.py"
        with open(fit_cfg_n, 'w') as fit_cfg_f:
            printFitCfg(fit_cfg_f, it, "/tmp", inListB[nFit], finListB[nFit], "Barrel", nFit)

        fitSrc_n = srcPath + f"/Fit/submit_EB_{nFit}_iter_{it}.sh"
        with open(fitSrc_n, 'w') as fitSrc_f:
            destination_s = f"{eosPath}/{dirname}/iter_{it}/{NameTag}Barrel_{nFit}_{calibMapName}"
            logpath = f"{pwd}/{dirname}/log/fitEpsilonPlot_EB_{nFit}_iter_{it}.log"
            tmpdir = "$TMPDIR" if (isOtherT2 and storageSite=="T2_BE_IIHE" and isCRAB) else "/tmp"
            printSubmitFitSrc(
                fitSrc_f,
                fit_cfg_n,
                f"{tmpdir}/{NameTag}Barrel_{nFit}_{calibMapName}",
                destination_s,
                pwd,
                logpath
            )

        subprocess.run(['chmod', '777', fitSrc_n], check=True)

    for nFit in range(nEE):
        fit_cfg_n = cfgFitPath + f"/fitEpsilonPlot_EE_{nFit}_iter_{it}.py"
        with open(fit_cfg_n, 'w') as fit_cfg_f:
            printFitCfg(fit_cfg_f, it, "/tmp", inListE[nFit], finListE[nFit], "Endcap", nFit)

        fitSrc_n = srcPath + f"/Fit/submit_EE_{nFit}_iter_{it}.sh"
        with open(fitSrc_n, 'w') as fitSrc_f:
            destination_s = f"{eosPath}/{dirname}/iter_{it}/{NameTag}Endcap_{nFit}_{calibMapName}"
            logpath = f"{pwd}/{dirname}/log/fitEpsilonPlot_EE_{nFit}_iter_{it}.log"
            tmpdir = "$TMPDIR" if (isOtherT2 and storageSite=="T2_BE_IIHE" and isCRAB) else "/tmp"
            printSubmitFitSrc(
                fitSrc_f,
                fit_cfg_n,
                f"{tmpdir}/{NameTag}Endcap_{nFit}_{calibMapName}",
                destination_s,
                pwd,
                logpath
            )

        subprocess.run(['chmod', '777', fitSrc_n], check=True)

# build command with options and arguments
calibCMD = f"python3 {pwd}/calibJobHandlerCondor.py -n {ijob}"
if options.recoverFill: calibCMD += " --recover-fill"
if options.daemonLocal: calibCMD += " --daemon-local"
if options.tokenFile:   calibCMD += f" --token-file {options.tokenFile}"
if options.minEfficiencyToRecoverFill >= 0.0:
    calibCMD += f" --min-efficiency-recover-fill {options.minEfficiencyToRecoverFill}"
calibCMD += "\n"

# generate submit script
env_script_n = f"{workdir}/submit.sh"
with open(env_script_n, 'w') as env_script_f:
    env_script_f.write("#!/bin/bash\n")
    env_script_f.write(f"cd {pwd}\n")
    env_script_f.write("ulimit -c 0\n")
    env_script_f.write("eval `scramv1 runtime -sh`\n")
    env_script_f.write(calibCMD)
    env_script_f.write(f"rm -rf {pwd}/core.*\n")

subprocess.run(['chmod', '777', env_script_n])

# dummy exec for condor
with open(f"{condordir}/dummy_exec_daemon.sh", 'w') as dummy_exec:
    dummy_exec.write('#!/bin/bash\n')
    dummy_exec.write('bash $*\n')

# condor submit file
condor_file_name = f"{condordir}/condor_submit_daemon.condor"
with open(condor_file_name, 'w') as condor_file:
    condor_file.write(f'''Universe = vanilla
Executable = {os.path.abspath(condordir)}/dummy_exec_daemon.sh
use_x509userproxy = True
x509userproxy = $ENV(X509_USER_PROXY)
Log        = {os.path.abspath(condordir)}/$(ProcId).log
Output     = {os.path.abspath(condordir)}/$(ProcId).out
Error      = {os.path.abspath(condordir)}/$(ProcId).error
getenv      = True
environment = "LS_SUBCWD={os.environ['PWD']}"
request_memory = 4000
+MaxRuntime = 604800
+JobBatchName = "ecalpro_daemon"
''')
    if os.environ['USER'] in ['mciprian']:
        condor_file.write('+AccountingGroup = "group_u_CMS.CAF.ALCA"\n\n')
    else:
        condor_file.write('\n')
    condor_file.write(f'arguments = {os.path.abspath(env_script_n)} \nqueue 1 \n\n')

# submit jobs
submit_s = f'condor_submit {condor_file_name}'
if not options.create:
    print(f"[calib] Number of jobs created = {ijob}")
    print("[calib] Submitting calibration handler")
    if options.daemonLocal:
        print(f"[calib]  '-- source {os.path.abspath(env_script_n)}'")
        os.system(f"source {os.path.abspath(env_script_n)}")
    else:
        print(f"[calib]  '-- {submit_s}'")
        output = subprocess.run(submit_s, shell=True, capture_output=True, text=True)
        print(f"[calib]  '-- {output.stdout.splitlines()[0]}'")
else:
    print("options -c was given: jobs are not submitted, but all folders and files were created normally. You can still do local tests.")
    print(f"To run the whole code use the following command:\n{submit_s}")

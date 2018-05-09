#!/usr/bin/python                                                                                                                                                             

import os
import time

currentDir = os.getcwd()
print ("currentDir = %s" %currentDir)


nRUN = 6
#print 'number of jobs = %d' %(len(nRUN))
print ('number of jobs = %d' %nRUN)

#inputCard = "delphes_card_HLLHC.tcl"
#inputCard = "delphes_card_HLLHC_v1.tcl"
#inputCard = "delphes_card_HLLHC_v2.tcl"
#inputCard = "delphes_card_HLLHC_v3.tcl"
#inputCard = "delphes_card_HLLHC_vall.tcl"
#inputCard = "delphes_card_HLLHC_vreso.tcl"
#inputCard = "delphes_card_HLLHC_veff.tcl"
#inputCard = "delphes_card_HLLHC_vd0.tcl"
#inputCard = "delphes_card_HLLHC_vdz.tcl"
#inputCard = "delphes_card_HLLHC_vd0Lot.tcl"
inputCard = "delphes_card_HLLHC_vd0x5.tcl"
#inputCard = "delphes_card_HLLHC_vd0x2.tcl"
#inputCard = "delphes_card_HLLHC_vd0x1p5.tcl"

#inputPhytia = "configNoLHE.cmnd"
#QCD
inputPhytia = "configNoLHE_qcd_pT160.cmnd"
#inputPhytia = "configNoLHE_qcd.cmnd"


outDir = "QCD_vd0x5_pT160"
#outDir = "QCD_vd0x5_pT140"
#outDir = "QCD_vd0x5_pT120"
#outDir = "QCD_vd0x5_pT100"
#outDir = "QCD_vd0x5_pT80"

#outDir = "QCD_v0_pT160"
#outDir = "QCD_v0_pT140"
#outDir = "QCD_v0_pT120"
#outDir = "QCD_v0_pT100"
#outDir = "QCD_v0_pT80"

#outDir = "QCD_vd0x2_pT160"
#outDir = "QCD_vd0x2_pT140"
#outDir = "QCD_vd0x2_pT120"
#outDir = "QCD_vd0x2_pT100"
#outDir = "QCD_vd0x2_pT80"



#outDir = "TTB_v0_b"
#outDir = "QCD_v0_b"

#outDir = "QCD_v0_highStat_13Mar"

#outDir = "TTB_vd0x5_b"
#outDir = "QCD_vd0x5_b"

eosArea = "/eos/cms/store/user/amartell/Pheno/"
#eosArea = "/eos/user/a/amartell/Delphes/"

subDir1 = currentDir+"/"+outDir;
os.system('mkdir %s' %(subDir1));

#subDir2 = "/tmp/amartell/"+outDir;
#os.system('mkdir %s' %(subDir2));

eosDir = eosArea+"/"+outDir;
os.system('eos  mkdir -p %s' %(eosDir));

outLancia = open('%s/lanciaProd_%s.sh' %(currentDir, outDir),'w');

for num in range(0+nRUN):
    print (num);
    subDir = currentDir+"/"+outDir+"/JOB_%d" %(num);
    os.system('mkdir %s' %(subDir));
    print ('subDir = %s  >>>> DID YOU CREATE THIS?' %subDir);

    outScript = open('%s/bjob.sh' %(subDir),'w');
    outScript.write('#!/bin/bash \n');
    outScript.write('cd %s \n' %(subDir));
    outScript.write('export SCRAM_ARCH=slc6_amd64_gcc530 \n');
    outScript.write('cd ../../../../CMSSW_9_1_0_pre3/  \n');
    outScript.write('eval `scramv1 ru -sh` \n');
    outScript.write('cd %s \n' %(currentDir));
    outScript.write('cd .. \n');
    #outScript.write('cmsMkdir %s \n' %(outFolder));
    outScript.write('pwd \n ')
    outROOTFileName = "/tmp/amartell/%s_delphes_HLLHCwithSmearing_PFJ_10k_%d.root" %(outDir, num);
    outScript.write('./DelphesPythia8 cards/%s examples/Pythia8/%s %s \n' %(inputCard, inputPhytia, outROOTFileName) );
    outScript.write('xrdcp -N -v %s root://eoscms.cern.ch/%s/' %(outROOTFileName, eosDir));
    #outScript.write('cp -N -v %s %s/%s' %(outROOTFileName, eosDir, outROOTFileName));
    outScript.write( '\n ' );
    os.system('chmod 777 %s/bjob.sh' %(subDir));
    outLancia.write(' bsub -cwd %s -q 8nh %s/bjob.sh \n' %(subDir, subDir));
    time.sleep(0.5); # seconds
    os.system(' bsub -cwd %s -q 8nh %s/bjob.sh \n' %(subDir, subDir));

# RGA_first_exp_dihadron

Description of the java code used by Timothy Hayward for the RGA first experiemnt dihadron analysis. Full analysis note and rough draft of the PRL letter are located in https://clas12-docdb.jlab.org/cgi-bin/DocDB/private/ShowDocument?docid=744. 

   hayward_coatjava_extensions.jar
   
Place this jar in "myClara/plugins/clas12/lib/services/" in order to use the java packages "extended_kinematic_fitters" and "analyzers".

   analysis_fitter.java
   
This is the extension to the CLAS12 eventbuilder that I use to identify electrons and pions. It first identifies electrons and pions assigned by the EventBuilder and then adds in an extra layer on top of that of so-called PID enhancements and fiducial cuts.

   HadronPairKinematics.java
   
This is the class I use to identify and return all relevant kinematic variables for the dihadron analysis. It includes the boolean "channel_test" that I use to determined whether an event passes all of our channel selection cuts described in the analysis note. 

   RGA_dihadrons_cross_check.groovy
   
Example groovy script used to create the cross check text file here. Shows how to access the relevant kinematic variables using the kinematic fitter and dihadron analysis code. I placed the skim4 run_005481 in "/work/clas12/thayward/temp/005418/" (although I do not intend to maintain this for long). You can run this script by indicating a directory with hipo files inside of it as the argument.

   run_005481_kinematic_comparisons.txt
   
Text file containing the output of the variables listed below for the skim4 run_005481. At the time of writing this was located in "/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass1/v0/dst/train/skim4/skim4_005418.hipo". Some event numbers may repeat but this is when there are 3 or more pions (and thus multiple pairs) in the event.

     Event#, e_p, pi+_p, pi-_p, Q2, W, MX, xF1, xF2, y, x, z, Mh, phperp, phiH, phiR, theta


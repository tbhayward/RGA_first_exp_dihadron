/*
 * author Timothy B. Hayward
 * 
 * RGA first publication dihadron cross check
 */

import java.io.File;

import org.jlab.io.hipo.*;
import org.jlab.io.base.DataEvent;
import org.jlab.clas.physics.*;
import org.jlab.clas12.physics.*;

// import from RGA_dihadrons
import extended_kinematic_fitters.*; 
import analyzers.*;


public class RGA_dihadrons_cross_check {

	public static void main(String[] args) {

		File[] hipo_list;
		if (args.length == 0) {
			// exits program if input directory not specified 
    	   	println("ERROR: Please enter a hipo file directory as the first argument");
    	  	System.exit(0);
    	} else {
    		File directory = new File(args[0]);
    		hipo_list = directory.listFiles();
    	}
    	int n_files;
		if ((args.length < 2)||(Integer.parseInt(args[1])>hipo_list.size())) {
			// if number of files not specified or too large, set to number of files in directory
			println("WARNING: Number of files not specified or number too large."); 
			println("Setting # of files to be equal to number of files in the directory.");
			n_files = hipo_list.size();
		} else{
			// if specified, convert to int
			n_files = Integer.parseInt(args[1]);
		}

		GenericKinematicFitter analysis_fitter = new analysis_fitter(10.6041);
		EventFilter filter = new EventFilter("11:211:-211:X+:X-:Xn");

		for (int current_file; current_file<n_files; current_file++) {
			println(); println(); println("Opening file "+Integer.toString(current_file+1)
				+" out of "+n_files);
			// limit to a certain number of files defined by n_files
			HipoDataSource reader = new HipoDataSource();
			reader.open(hipo_list[current_file]); // open next hipo file
			HipoDataEvent event = reader.getNextEvent(); 

			while(reader.hasEvent()==true){ // cycle through events
				event = reader.getNextEvent();

				PhysicsEvent research_Event = analysis_fitter.getPhysicsEvent(event);
				if(filter.isValid(research_Event)==true){
					int num_piplus = research_Event.countByPid(211); 
					int num_piminus = research_Event.countByPid(-211);

					// cycle over multiple hadron pairs if they exist
					for (int current_plus = 0; current_plus < num_piplus; current_plus++) {
						for (int current_minus = 0; current_minus < num_piminus; current_minus++) {

							HadronPairKinematics dihadronVars = new HadronPairKinematics(
								event, research_Event, current_plus, current_minus, false);

							if (dihadronVars.channel_test(dihadronVars)) {
								int evnum = event.getBank("RUN::config").getInt('event',0)
								int helicity = dihadronVars.get_helicity();

								// SIDIS cut variables
								double Q2 = dihadronVars.Q2();
								double W = dihadronVars.W();

								// exclusive region cut variables
								double MX = dihadronVars.missingMass();

								// current fragmentation region cut variables
								double xF1 = dihadronVars.xF1();
								double xF2 = dihadronVars.xF1();

								// radiative effect cut variables
								double y = dihadronVars.y();

								// BSA variables
								double x = dihadronVars.x();
								double z = dihadronVars.z();
								double mpipi = dihadronVars.pairMass();								
								double pHT = dihadronVars.pHT(); // phperp

								// angles
								double phiH = dihadronVars.phiH();
								double phiR = dihadronVars.phiR();
								double theta = dihadronVars.theta();

								System.out.print(evnum)
								System.out.printf(", %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, ",
									Q2, W, MX, xF1, xF2, y);
								System.out.printf("%.2f, %.2f, %.2f, %.2f, ",
									x, z, mpipi, pHT);
								System.out.printf("%.2f, %.2f, %.2f",
									phiH, phiR, theta);
								println();
							}
						}
					}
				}
			}
		}
	}
}
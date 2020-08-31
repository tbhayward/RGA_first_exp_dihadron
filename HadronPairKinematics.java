package analyzers;
/**
 *
 * by timothy b. hayward, with anselm vossen
 */

import org.jlab.clas.physics.Particle;
import org.jlab.clas.physics.PhysicsEvent;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataBank;

import org.jlab.clas.physics.*;

public class HadronPairKinematics {
    
    protected byte helicity;
    
    protected int num_elec;
    protected int num_piplus;
    protected int num_piminus;
    protected int num_kplus;
    protected int num_kminus;
    protected int num_particles;
    
    protected double Q2;
    protected double nu;
    protected double x;
    protected double W;
    protected double y;   
    protected double z;
    protected double z1;
    protected double z2;
    protected double gamma;
    
    protected double Depolarization_A;
    protected double Depolarization_B;
    protected double Depolarization_C;
    protected double Depolarization_V;
    protected double Depolarization_W;
    
    protected double pairMass;
    protected double pT;
    protected double pHT;
    protected double xF;
    protected double xF1;
    protected double xF2;
    protected double theta;
    protected double phiH;
    protected double phiR;
    
    protected double elec_p;
    protected double piplus_p, piminus_p, pipi_p, boosted_pipCOM_p, boosted_pipi_p, dot_product;
    protected double elec_theta;
    protected double missingMass;
    
    public int get_helicity() {
        return helicity; // return helicity of event, 0 or 1
    }
    
    public int num_elec() {
        return num_elec; // returns number of electrons
    }
    
    public int num_piplus() {
        return num_piplus; // returns number of piplus
    }
    
    public int num_piminus() {
        return num_piminus; // returns number of piminus
    }
    
    public int num_kplus() {
        return num_kplus; // returns number of kplus
    }
    
    public int num_kminus() {
        return num_kminus; // returns number of kminus
    }
    
    public double Q2() {
        return Q2; // returns Q2
    }
    
    public double nu() {
        return nu; // returns nu
    }
    
    public double x() {
        return x; // returns x
    }
    
    public double W() {
        return W; // returns W
    }
    
    public double y() {
        return y; // returns y
    }
    
    public double z() {
        return z; // returns z
    }
    
    public double z1() {
        return z1; // returns z1
    }
    
    public double z2() {
        return z2; // returns z2
    }
    
    public double gamma() {
        return gamma; // returns z2
    }
    
    public double Depolarization_A() {
        return Depolarization_A; // 
    }
    
    public double Depolarization_B() {
        return Depolarization_B; // 
    }
    
    public double Depolarization_C() {
        return Depolarization_C; // 
    }
    
    public double Depolarization_V() {
        return Depolarization_V; // 
    }
    
    public double Depolarization_W() {
        return Depolarization_W; // 
    }
    
    public double pairMass() {
        return pairMass;
    }
    
    public double pT() {
        return pT;
    }
    
    public double pHT() {
        return pHT;
    }
    
    public double xF() {
        return xF;
    }
    
    public double xF1() {
        return xF1;
    }
    
    public double xF2() {
        return xF2;
    }
    
    public double theta() {
        return theta;
    }
    
    public double phiH() {
        return phiH;
    }
    
    public double phiR() {
        return phiR;
    }
    
    public double elec_p() {
        return elec_p;
    }
    
    public double elec_theta() {
        return elec_theta;
    }
    
    public double piplus_p() {
        return piplus_p;
    }
    
    public double piminus_p() {
        return piminus_p;
    }
    
    public double pipi_p() {
        return pipi_p;
    }
    
    public double boosted_pipCOM_p() {
        return boosted_pipCOM_p;
    }
    
    public double boosted_pipi_p() {
        return boosted_pipi_p;
    }
    
    public double dot_product() {
        return dot_product;
    }
    
    public double missingMass() {
        return missingMass;
    }
    
    
    private static double particle_mass (int pid) {
	if (pid==11) {
            return 0.0005109989461;
	} else if (pid==211||pid==-211) {
            return 0.139570;
	} else if (pid==321||pid==-321) {
            return 0.493677;
	} else if (pid==2212||pid==-2212) {
            return 0.938272;
	}
            return -1;
    }
    
    public static boolean channel_test(HadronPairKinematics dihadronVars) {
	boolean channel_success = true;
//        if (dihadronVars.helicity==0){
//            channel_success = false;
//        }
	if (dihadronVars.Q2()<1) {
            channel_success = false;
	} else if (dihadronVars.W()<2) {
            channel_success = false;
	} 
//////        else if (dihadronVars.z1()<0.0 || dihadronVars.z2()<0.0) {
//////            channel_success = false;
// 
        else if (dihadronVars.xF1()<0.0 || dihadronVars.xF2()<0.0) {
            channel_success = false;
	} 
        else if (dihadronVars.y()>0.8) {
            channel_success = false;
        } 
//////        else if (dihadronVars.z()>0.95 || dihadronVars.z()<0.1) {
//////            channel_success = false;
//////	} 
        else if (dihadronVars.missingMass()<1.5) {
            channel_success = false;
	} 
//        else if (dihadronVars.missingMass()<1.5 && dihadronVars.missingMass()>1.15) {
//            channel_success = false;
//	} 
//        else if (dihadronVars.pairMass()<0.63) {
//            channel_success = false;
//        }
//        else if (dihadronVars.pairMass()<0.63 || dihadronVars.pairMass()>0.89) {
//            channel_success = false;
//        }
//        else if (dihadronVars.pairMass()<0.89) {
//            channel_success = false;
//        }
//        else if (dihadronVars.Q2() > 2.14) {
//            channel_success = false;
//        }
//        else if (dihadronVars.Q2() < 2.140 || dihadronVars.Q2() > 3.10) {
//            channel_success = false;
//        }
//        else if (dihadronVars.Q2() < 3.10 ) {
//            channel_success = false;
//        }
//        else if (dihadronVars.pHT() > 0.360) {
//            channel_success = false;
//        }
//        else if (dihadronVars.pHT() < 0.360 || dihadronVars.pHT() > 0.585) {
//            channel_success = false;
//        }
//        else if (dihadronVars.pHT() < 0.585 ) {
//            channel_success = false;
//        }
	return channel_success;
	}
    
    public HadronPairKinematics(DataEvent event, PhysicsEvent recEvent, int plusIndex, int minusIndex, boolean isMC) {
        
        if (!isMC) {
            HipoDataBank eventBank = (HipoDataBank) event.getBank("REC::Event");
            HipoDataBank recBank = (HipoDataBank) event.getBank("REC::Particle"); 
        
            helicity = eventBank.getByte("helicity", 0);
        }
        
	num_elec = recEvent.countByPid(11); // returns number of electrons
	num_piplus = recEvent.countByPid(211); 
	num_piminus = recEvent.countByPid(-211);
	num_kplus = recEvent.countByPid(321);
	num_kminus = recEvent.countByPid(-321);
        num_particles = num_elec+num_piplus+num_piminus+num_kplus+num_kminus;
        
        LorentzVector lv_beam = new LorentzVector();
	lv_beam.setPxPyPzM(0, 0, Math.pow(10.6041*10.6041-particle_mass(11)*particle_mass(11),0.5), 
                particle_mass(11));
	LorentzVector lv_target = new LorentzVector();
	lv_target.setPxPyPzM(0, 0, 0, particle_mass(2212));
        String electron_index = "[11,0]"; // highest p 
	Particle scattered_electron= recEvent.getParticle(electron_index); // 
        elec_p = scattered_electron.p();
        elec_theta = scattered_electron.theta();
	LorentzVector lv_e = new LorentzVector();
	lv_e.setPxPyPzM(scattered_electron.px(), scattered_electron.py(), 
            scattered_electron.pz(), particle_mass(11));
        
        LorentzVector q = new LorentzVector(lv_beam); q.sub(lv_e);
	Q2 = -q.mass2();
	nu = lv_beam.e()-lv_e.e();
	x  = Q2 / (2 * particle_mass(2212) * nu);
	W  = Math.pow(Math.pow(particle_mass(2212),2)+2*particle_mass(2212)*nu - Q2, 0.5);
	y = nu/lv_beam.e();
        double M = 0.938272088; // mass of proton
        gamma = 2*M*x/Math.pow(Q2, 0.5);
        
        Depolarization_A = 1/(1+gamma*gamma)*(1-y+y*y/2+y*y*gamma*gamma/4);
        Depolarization_B = 1/(1+gamma*gamma)*(1-y-y*y*gamma*gamma/4);
        Depolarization_C = (y/Math.pow(1+gamma*gamma, 0.5))*(1-y/2);
        Depolarization_V = (2-y)/(1+gamma*gamma)*Math.pow(1-y-y*y*gamma*gamma/4,0.5);
        Depolarization_W = y/(Math.pow(1+gamma*gamma, 0.5))*Math.pow(1-y-y*y*gamma*gamma/4,0.5);
        
        LorentzVector gN = new LorentzVector(q);
	gN.add(lv_target);
	Vector3 gNBoost = gN.boostVector();
	gNBoost.negative();
       
        DISKinematics(gN, gNBoost, recEvent, q, lv_e, lv_target, plusIndex, minusIndex);
    }
    
    public final void DISKinematics(LorentzVector gN, Vector3 gNBoost, PhysicsEvent recEvent, 
            LorentzVector q, LorentzVector lv_e, LorentzVector lv_target, int plusIndex, int minusIndex) {
        String piplus_index = "[211,"+plusIndex+"]";
        String piminus_index = "[-211,"+minusIndex+"]";
        String combined_index = "[211,"+plusIndex+"]+[-211,"+minusIndex+"]";
        Particle pipi_particle = recEvent.getParticle(combined_index);
        pairMass = pipi_particle.mass();
        
        LorentzVector pipi = new LorentzVector();
        pipi.setPxPyPzM(pipi_particle.px(), pipi_particle.py(), pipi_particle.pz(), pipi_particle.mass());
        Particle piplus_particle = recEvent.getParticle(piplus_index);
        pipi_p = pipi.p();
        
        LorentzVector piplus = new LorentzVector();
	piplus.setPxPyPzM(piplus_particle.px(), piplus_particle.py(), piplus_particle.pz(), piplus_particle.mass());
        piplus_p = piplus.p();
	Particle piminus_particle = recEvent.getParticle(piminus_index);
	LorentzVector piminus = new LorentzVector();
	piminus.setPxPyPzM(piminus_particle.px(), piminus_particle.py(), piminus_particle.pz(), 
	piminus_particle.mass());
        piminus_p = piminus.p();
        
        LorentzVector mx = new LorentzVector(q); mx.add(lv_target); mx.sub(piplus); mx.sub(piminus);
        missingMass = mx.mass(); 
        
        // boost to Breit frame
	LorentzVector boosted_pipi = new LorentzVector(pipi);
	boosted_pipi.boost(gNBoost);
	LorentzVector boosted_piplus = new LorentzVector(piplus);
	boosted_piplus.boost(gNBoost);
	LorentzVector boosted_piminus = new LorentzVector(piminus);
	boosted_piminus.boost(gNBoost);
	LorentzVector boosted_q = new LorentzVector(q);
	boosted_q.boost(gNBoost); 
	Vector3 boosted_q_unit = new Vector3(); 
	boosted_q_unit.setMagThetaPhi(1, boosted_q.theta(), boosted_q.phi());
	LorentzVector boosted_lv_e = new LorentzVector(lv_e);
	boosted_lv_e.boost(gNBoost); 
	Vector3 boosted_lv_e_unit = new Vector3(); 
	boosted_lv_e_unit.setMagThetaPhi(1, boosted_lv_e.theta(), boosted_lv_e.phi());
        
        LorentzVector piCOM = new LorentzVector(pipi);
        Vector3 piCOMBoost = piCOM.boostVector();
        piCOMBoost.negative();
        LorentzVector boosted_pipCOM = new LorentzVector(piplus); boosted_pipCOM.boost(piCOMBoost);
        LorentzVector boosted_pimCOM = new LorentzVector(piminus); boosted_pimCOM.boost(piCOMBoost);
        LorentzVector boosted_pipiCOM = new LorentzVector(pipi); boosted_pipiCOM.boost(piCOMBoost);
        boosted_pipCOM_p = boosted_pipCOM.p();
        boosted_pipi_p = boosted_pipi.p();

        dot_product = boosted_pipCOM.vect().dot(pipi.vect());
        double costheta = (boosted_pipCOM.vect().dot(pipi.vect())) / 
                (boosted_pipCOM.vect().mag()*pipi.vect().mag());
        theta = Math.acos(costheta);

        pT = boosted_q_unit.cross(boosted_pipi.vect()).mag();
	pHT = boosted_q_unit.cross(boosted_pipi.vect()).mag();
//        pHT = boosted_q_unit.dot(boosted_pipi.vect());
        
        
        xF = 2*(boosted_pipi.vect().dot(boosted_q.vect())) /(boosted_q.vect().mag()*W);
        xF1 = 2*(boosted_piplus.vect().dot(boosted_q.vect())) /(boosted_q.vect().mag()*W);
        xF2 = 2*(boosted_piminus.vect().dot(boosted_q.vect())) /(boosted_q.vect().mag()*W);
        
        z = pipi.e()/q.e();
        
        z1 = piplus.e()/q.e();
	z2 = piminus.e()/q.e();

	Vector3 pipiBoostVect = boosted_pipi.boostVector();
	pipiBoostVect.negative();
	LorentzVector boosted_boosted_piplus = new LorentzVector(piplus.px(), piplus.py(), piplus.pz(), 
	piplus.e());
	boosted_boosted_piplus.boost(pipiBoostVect);
//	theta = Math.acos(boosted_boosted_piplus.vect().dot(boosted_pipi.vect()) / ( 
//            boosted_boosted_piplus.vect().mag()* boosted_pipi.vect().mag()));
        
	Vector3 vecR = new Vector3(boosted_piplus.vect()); // not really R yet
	vecR.setMagThetaPhi(boosted_piplus.vect().mag()/z1, boosted_piplus.vect().theta(), boosted_piplus.vect().phi());
	Vector3 vecH = new Vector3();
	vecH.setMagThetaPhi(boosted_piminus.vect().mag()/z2, boosted_piminus.vect().theta(),
            boosted_piminus.vect().phi());
	vecR.sub(vecH); // this is really R now that we do this subtraction

	Vector3 vectRt = new Vector3();
	Vector3 R_Q = new Vector3();
							
	R_Q.setMagThetaPhi(vecR.dot(boosted_q_unit), boosted_q_unit.theta(), boosted_q_unit.phi());
	vectRt=vecR;
	vectRt.sub(R_Q);

	Vector3 vectPh = new Vector3(boosted_pipi.vect());
	Vector3 Pt_Q = new Vector3();
	Pt_Q.setMagThetaPhi(vecR.dot(boosted_q_unit), boosted_q_unit.theta(), boosted_q_unit.phi());
	Vector3 vectPhT = new Vector3(vectPh);
	vectPhT.sub(Pt_Q);
//        pHT = vectPhT.mag();
							
	Vector3 vT = new Vector3(boosted_q_unit.cross(boosted_lv_e.vect())); vT.unit();
	Vector3 vTR = new Vector3(boosted_q_unit.cross(vectRt)); vTR.unit();
	Vector3 vTH = new Vector3(boosted_q_unit.cross(vectPhT)); vTH.unit();

	double cosPhiR = vT.dot(vTR);
	double sinPhiR = boosted_lv_e.vect().cross(vectRt).dot(boosted_q_unit);
	double cosPhiH = vT.dot(vTH);
	double sinPhiH = boosted_lv_e.vect().cross(vectPhT).dot(boosted_q_unit);

	// scaling
	double rScale = boosted_q_unit.cross(boosted_lv_e.vect()).mag()*boosted_q_unit.cross(vectRt).mag();
	sinPhiR = sinPhiR/rScale;
	double hScale = boosted_q_unit.cross(boosted_lv_e.vect()).mag()*boosted_q_unit.cross(vectPh).mag();
	sinPhiH = sinPhiH/hScale;
	phiH = Math.acos(cosPhiH);
	phiR = Math.acos(cosPhiR);
	if (sinPhiR < 0.0) {
            phiR = 2*Math.PI - phiR;
	}
	if (sinPhiH < 0.0) {
            phiH = 2*Math.PI - phiH;
	}
               
    }
}
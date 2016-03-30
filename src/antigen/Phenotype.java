package antigen;
/* Antigenic phenotype present in individual Viruses and within Hosts as immune history */
/* Should be able to calculate distance and cross-immunity between two phenotypes */
/* Moving up to multiple dimensions is non-trivial and requires thought on the implementation */
/* Multiple Viruses can reference a single Phenotype object */

import static java.lang.Math.*;
import java.util.*;

public class Phenotype {

	// fields
	private float traitA;
	private float traitB;	
	
	// constructor
	public Phenotype() {
	
	}
	public Phenotype(double tA, double tB) {
		traitA = (float)tA;
		traitB = (float)tB;
	}
	
	// mean constructor
	public Phenotype(List<Phenotype> samples) {
		double traitASum = 0;
		double traitBSum = 0;
		for(Phenotype sample : samples) {
			traitASum += sample.traitA;
			traitBSum += sample.traitB;
		}
		traitA = (float)(traitASum / samples.size());
		traitB = (float)(traitBSum / samples.size());
	}
	
	public Phenotype mutateNormal(double sd)
	{
		// Isotropic normal error
		double theta = Random.nextDouble(0,2*Math.PI);
		double r = Random.nextNormal(0.0, sd);
		return new Phenotype(traitA + r * Math.cos(theta), traitB + r * Math.sin(theta));
	}
	
	public double getTraitA() {
		return traitA;
	}
	public double getTraitB() {
		return traitB;
	}	
	
	public void setTraitA(double tA) {
		traitA = (float)tA;
	}
	public void setTraitB(double tB) {
		traitB = (float)tB;
	}		
		
	// raw antigenic distance between two phenotypes
	public double distance(Phenotype p) {
		double distA = (getTraitA() - p.getTraitA());
		double distB = (getTraitB() - p.getTraitB());	
		double dist = (distA * distA) + (distB * distB);
		dist = Math.sqrt(dist);
		return dist;
	}

	// cross immunity between a virus phenotype and a host's immune history
	// here encoded more directly as risk of infection, which ranges from 0 to 1
	public double riskOfInfection(
			List<Phenotype> immuneHistory, List<Phenotype> vaccinationHistory,
			double smithConversion, double homologousImmunity, double vaccineImmuneBreadth
	) {
	
		// find closest phenotype in history
		double closestDistance = 100.0;
		for(int i = 0; i < immuneHistory.size(); i++) {
			double thisDistance = distance(immuneHistory.get(i));
			if(thisDistance < closestDistance) {
				closestDistance = thisDistance;
			}
//			if (thisDistance < 0.01) {
//				break;
//			}
		}
//		if (immuneHistory.size() > 0) {
//			for (int i = immuneHistory.size() - 1; i >= 0; i--) {
//				double thisDistance = distance(immuneHistory.get(i));
//				if (thisDistance < closestDistance) {
//					closestDistance = thisDistance;
//				}
//			}
//		}
		
		double risk = (closestDistance == Double.POSITIVE_INFINITY) ?
				1.0 : (closestDistance * smithConversion);
		double minRisk = 1.0 - homologousImmunity;
		risk = Math.max(minRisk, risk);
		risk = Math.min(1.0, risk);
		
		if(vaccinationHistory != null) {
			double closestVacDistance = Double.POSITIVE_INFINITY;
			for(int i = vaccinationHistory.size() - 1; i >= 0; i--) {
				Phenotype vaccineStrain = vaccinationHistory.get(i);
				double distance = distance(vaccineStrain);
				if(distance < closestVacDistance) {
					closestVacDistance = distance;
				}
			}
			if(closestVacDistance != Double.POSITIVE_INFINITY) {
				double vaccineRisk = closestVacDistance * smithConversion / vaccineImmuneBreadth;
				vaccineRisk = Math.max(minRisk, vaccineRisk);
				vaccineRisk = Math.min(1.0, vaccineRisk);
				
				risk = min(risk, vaccineRisk);
			}
		}
		
		return risk;
		
	}
		
	// returns a mutated copy, original Phenotype is unharmed
	// mutate with gamma
	public Phenotype mutate(boolean mut2D, double meanStep, double sdStep) {
		
		// direction of mutation
		double theta = 0;
		if (mut2D) {
			theta = Random.nextDouble(0,2*Math.PI);
		} else {
			if (Random.nextBoolean(0.5)) { theta = 0; }
			else { theta = Math.PI; }
		}
		
		// size of mutation
		double r = meanStep;
		if (sdStep != 0.0) {
			double alpha = (meanStep *  meanStep) / (sdStep * sdStep);
			double beta = (sdStep * sdStep) / meanStep;
			r = Random.nextGamma(alpha, beta);
		}
		
		// create phenotype
		double mutA = getTraitA() + r * Math.cos(theta);
		double mutB = getTraitB() + r * Math.sin(theta);
		Phenotype mutP = new Phenotype(mutA,mutB);
		return mutP;
				
	}	
	
	public String toString() {
		String fullString = String.format("%.4f,%.4f", traitA, traitB);
		return fullString;
	}
	
	public boolean equals(Object o) {
		if(!(o instanceof Phenotype)) {
			return false;
		}
		Phenotype gp = (Phenotype)o;
		return gp.traitA == traitA && gp.traitB == traitB;
	}
	
}

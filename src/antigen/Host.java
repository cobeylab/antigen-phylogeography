package antigen;
/* A human individual that harbors viruses and immunity */

import java.util.*;
import java.io.*;

public class Host {

	// fields
	private Virus infection;												
	private List<Phenotype> immuneHistory = new ArrayList<Phenotype>(0);
	private IterableBitSet vaccinationHistory;
	
	// Naive host
	public Host(boolean vaccinationActive) {
		if(vaccinationActive) {
			vaccinationHistory = new IterableBitSet(0);
		}
	}
	
	// Immune host
	public Host(Phenotype immunity, boolean vaccinationActive) {
		assert immunity != null;
		addToImmuneHistory(immunity);
		if(vaccinationActive) {
			vaccinationHistory = new IterableBitSet(0);
		}
	}
	
	// Infected host
	public Host(Virus v, boolean vaccinationActive) {
		infection = v;
		if(vaccinationActive) {
			vaccinationHistory = new IterableBitSet(0);
		}
	}
	
	public void addToImmuneHistory(Phenotype p) {
		for(Phenotype pi : immuneHistory) {
			if(pi.equals(p)) {
				return;
			}
		}
		immuneHistory.add(p);
	}
	
	// infection methods
	public void reset() {
		infection = null;
		immuneHistory = new ArrayList<Phenotype>(0);
		if(vaccinationHistory != null) {
			vaccinationHistory = new IterableBitSet(0);
		}
	}
	
	public boolean isInfected() {
		boolean infected = false;
		if (infection != null) {
			infected = true;
		}
		return infected;
	}
	public Virus getInfection() {
		return infection;
	}
	public void infect(double time, Virus pV, int d) {
		Virus nV = new Virus(time, pV, d);
		infection = nV;
	}
	public void clearInfection() {
		Phenotype p = infection.getPhenotype();
		addToImmuneHistory(p);
		infection = null;
	}
	public int getHistoryLength() {
		return immuneHistory.size();
	}
	
	// make a new virus with the mutated phenotype
	public void mutate(boolean mut2D, double meanStep, double sdStep, double birth) {
		Virus mutV = infection.mutate(mut2D, meanStep, sdStep, birth);
		infection = mutV;
	}
	
	public void vaccinate(int vaccineId) {
		vaccinationHistory.set(vaccineId);
	}
	
	// history methods
	
	public List<Phenotype> getHistory() {
		return immuneHistory;
	}	
	
	public IterableBitSet getVaccinationHistoryIdSet(Simulation sim) {
		return vaccinationHistory;
	}
	
	public List<Phenotype> getVaccinationHistory(Simulation sim) {
		if(vaccinationHistory == null) {
			return null;
		}
		
		List<Phenotype> vhList = new ArrayList<Phenotype>(vaccinationHistory.cardinality());
		for(int vaccineId : vaccinationHistory) {
			vhList.add(sim.getVaccine(vaccineId));
		}
		return vhList;
	}
	
	public void printImmuneHistory() {
		for (int i = 0; i < immuneHistory.size(); i++) {
			System.out.println(immuneHistory.get(i));
		}
	}

	public void printInfection(PrintStream stream) {
		if (infection != null) {
			stream.print(infection.getPhenotype());
		}
		else {
			stream.print("n");
		}
	}
	
	public void printHistory(PrintStream stream) {
		if (immuneHistory.size() > 0) {
			stream.print(immuneHistory.get(0));
			for (int i = 1; i < immuneHistory.size(); i++) {
				stream.print(";" + immuneHistory.get(i));
			}
		}
		else {
			stream.print("n");
		}
	}	
		
	public String toString() {
		return Integer.toHexString(this.hashCode());
	}	
	
}
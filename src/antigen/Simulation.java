package antigen;
/* Simulation functions, holds the host population */

import java.util.*;
import java.io.*;

// import org.tmatesoft.sqljet.core.SqlJetException;
// import org.tmatesoft.sqljet.core.SqlJetTransactionMode;
// import org.tmatesoft.sqljet.core.table.*;

import com.google.gson.GsonBuilder;
// import com.javamex.classmexer.*;

public class Simulation {
	private Parameters params;
	
	// Time tracking
	int timestep;
	int endstep;
	double day;
	double vaccineSeasonEndDay;
	
	int updateVaccineEvery;
	int deployVaccineStart;
	boolean vaccineSeason;
	// Initial virus/immunity
	Virus urVirus = null;
	Phenotype urImmunity = null;	
	
	
	// fields
	private List<HostPopulation> demes = new ArrayList<HostPopulation>();
	private double diversity;
	private double tmrca;
	private double netau;
	private double serialInterval;
	private double antigenicDiversity;	
	
	private List<Double> diversityList = new ArrayList<Double>();
	private List<Double> tmrcaList = new ArrayList<Double>();	
	private List<Double> netauList = new ArrayList<Double>();		
	private List<Double> serialIntervalList = new ArrayList<Double>();		
	private List<Double> antigenicDiversityList = new ArrayList<Double>();
	private List<Double> nList = new ArrayList<Double>();
	private List<Double> sList = new ArrayList<Double>();	
	private List<Double> iList = new ArrayList<Double>();	
	private List<Double> tList = new ArrayList<Double>();		
	private List<Double> casesList = new ArrayList<Double>();			
	private List<Double> unexposedList = new ArrayList<Double>();
	
	private List<Phenotype> vaccineStrains;
	private double lastVaccineDate;
	private int deployedVaccine;
	private Phenotype idealVaccineStrain;
	
	private PrintStream vaccineStream;
	
// 	private SqlJetDb sampleDb;
	
	// constructor
	public Simulation(Parameters params) {
		this.params = params;
		
		if(params.randomSeed == null) {
			params.randomSeed = Random.seed();
		}
		else {
			Random.seed(params.randomSeed);
		}
		System.err.printf("seed: %d\n", params.randomSeed);
		System.out.printf("seed: %d\n", params.randomSeed);
		
		// Initialize time
		timestep = 0;
		endstep = (int) ((double) params.endDay/params.deltaT);
		day = timestep*params.deltaT;
		
		updateVaccineEvery = (int)(params.vaccineStep / params.deltaT);
		vaccineStrains = new ArrayList<Phenotype>();
		
		// Initialize urVirus/urImmunity
		urVirus = new Virus(getDate());
		urImmunity = new Phenotype(params.initialTraitA, 0);
		
		VirusTree.root = urVirus;
		VirusTree.params = params;
		VirusTree.sim = this;
		
		for (int i = 0; i < params.demeCount; i++) {
			HostPopulation hp = new HostPopulation(this, params, urImmunity, urVirus, i);
			demes.add(hp);
		}
		assert getI() > 0;
		
		params.initialIs = new int[params.demeCount];
		params.initialRs = new int[params.demeCount];
		for(int i = 0; i < params.demeCount; i++) {
			params.initialIs[i] = demes.get(i).initialI;
			params.initialRs[i] = demes.get(i).initialR;
		}
		
		// Write out read-in parameters
		try {
			PrintWriter jsonWriter = new PrintWriter("parameters_out.json");
			new GsonBuilder().setPrettyPrinting().create().toJson(params, jsonWriter);
			jsonWriter.close();
		}
		catch(Exception e) {
			throw new RuntimeException(e);
		}
		initSampleDb();
	}
	
	private void updateDeployedVaccine() {
		deployedVaccine = nVaccines()-1;
	}
	
	private void updateVaccineStrain() {
		List<Phenotype> samples = new ArrayList<>();
		for(HostPopulation deme : demes) {
			for(Host host : deme.infecteds) {
				Virus virus = host.getInfection();
				samples.add(virus.getPhenotype());
			}
		}
		lastVaccineDate = getDate();
		idealVaccineStrain = new Phenotype(samples);
		vaccineStrains.add(idealVaccineStrain.mutateNormal(params.vaccineSD));
		System.err.printf("vaccine strain: %s\n", vaccineStrains.get(vaccineStrains.size() - 1));
		printVaccine();
	}
	
	public int nVaccines() {
		return vaccineStrains.size();
	}
	
	public Phenotype getVaccine(int vaccineId) {
		return vaccineStrains.get(vaccineId);
	}
	
	public Phenotype getLastVaccineStrain() {
		return vaccineStrains.get(vaccineStrains.size() - 1);
	}
	
	// methods
	
	public int getUnexposed(){
		int count = 0;
		//TODO: reimplement if needed
//		for(int i = 0; i < params.demeCount; i++){
//			HostPopulation hp = demes.get(i);
//			count += hp.getUnexposed();
//		}
		return count;
	}
	
	public int getN() {
		int count = 0;
		for (int i = 0; i < params.demeCount; i++) {
			HostPopulation hp = demes.get(i);
			count += hp.getN();
		}
		return count;
	}
	
	public int getS() {
		int count = 0;
		for (int i = 0; i < params.demeCount; i++) {
			HostPopulation hp = demes.get(i);
			count += hp.getS();
		}
		return count;
	}	
	
	public int getI() {
		int count = 0;
		for (int i = 0; i < params.demeCount; i++) {
			HostPopulation hp = demes.get(i);
			count += hp.getI();
		}
		return count;
	}		
	
	public int getT() {
		int count = 0;
		for (int i = 0; i < params.demeCount; i++) {
			HostPopulation hp = demes.get(i);
			count += hp.getT();
		}
		return count;
 	}	
	
	public int getCases() {
		int count = 0;
		for (int i = 0; i < params.demeCount; i++) {
			HostPopulation hp = demes.get(i);
			count += hp.getCases();
		}
		return count;
	}	
	
	public double getDiversity() {
		return diversity;
	}		
	
	public double getNetau() {
		return netau;
	}	
	
	public double getTmrca() {
		return tmrca;
	}	
	
	public double getSerialInterval() {
		return serialInterval;	
	}		
	
	public double getAntigenicDiversity() {
		return antigenicDiversity;
	}		
		
	// proportional to infecteds in each deme
	public int getRandomDeme() {
		int n = Random.nextInt(0,getN()-1);
		int d = 0;
		int target = (demes.get(0)).getN();
		while (n < target) {
			d += 1;
			target += (demes.get(d)).getN();
		}
		return d;
	}
	
	// return random virus proportional to worldwide prevalence
	public Virus getRandomInfection() {
	
		Virus v = null;
		
		if (getI() > 0) {
	
			// get deme proportional to prevalence
			int n = Random.nextInt(0,getI()-1);
			int d = 0;
			int target = (demes.get(0)).getI();
			while (d < params.demeCount) {
				if (n < target) {
					break;
				} else {
					d++;
					target += (demes.get(d)).getI();
				}	
			}
			HostPopulation hp = demes.get(d);
					
			// return random infection from this deme
			if (hp.getI()>0) {
				Host h = hp.getRandomHostI();
				v = h.getInfection();
			}
		
		}
		
		return v;
		
	}
	
	// return random host from random deme
	public Host getRandomHost() {
		int d = Random.nextInt(0,params.demeCount-1);
		HostPopulation hp = demes.get(d);
		return hp.getRandomHost();
	}
	
	public double getAverageRisk(Phenotype p) {
		
		double averageRisk = 0;
		for (int i = 0; i < 10000; i++) {
			Host h = getRandomHost();
			List<Phenotype> immuneHistory = h.getHistory();
			List<Phenotype> vaccinationHistory = h.getVaccinationHistory(this);
			averageRisk += p.riskOfInfection(immuneHistory, vaccinationHistory,
					params.smithConversion, params.homologousImmunity, params.vaccineImmuneBreadth);
		}
		averageRisk /= 10000.0;
		return averageRisk;
		
	}
		
	public void printImmunity() {
	
		try {
			File immunityFile = new File("out.immunity");
			immunityFile.delete();
			immunityFile.createNewFile();
			PrintStream immunityStream = new PrintStream(immunityFile);
			
			for (double x = VirusTree.xMin; x <= VirusTree.xMax; x += 0.5) {
				for (double y = VirusTree.yMin; y <= VirusTree.yMax; y += 0.5) {
				
					Phenotype p = new Phenotype(x,y);
					double risk = getAverageRisk(p);
					immunityStream.printf("%.4f,", risk);
				
				}
				immunityStream.println();
			}
			
			immunityStream.close();
		} catch(IOException ex) {
			System.out.println("Could not write to file"); 
			System.exit(0);
		}
	
	}
	
	public void printHostPopulation() {
	
		try {
			File hostFile = new File("out.hosts");
			hostFile.delete();
			hostFile.createNewFile();
			PrintStream hostStream = new PrintStream(hostFile);
			for (int i = 0; i < params.demeCount; i++) {
				HostPopulation hp = demes.get(i);
				hp.printHostPopulation(hostStream);
			}
			hostStream.close();
		} catch(IOException ex) {
			System.out.println("Could not write to file"); 
			System.exit(0);
		}
	
	}	
		
	public void makeTrunk() {
		for (int i = 0; i < params.demeCount; i++) {
			HostPopulation hp = demes.get(i);
			hp.makeTrunk();
		}
	}
	
	public void printState() {
	
		System.out.printf("%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%d\t%d\t%d\t%d\t%d\n", (int) day, getDiversity(), getTmrca(),  getNetau(), getSerialInterval(), getAntigenicDiversity(), getN(), getS(), getI(), getT(), getCases(), getUnexposed());

// 		if (params.memoryProfiling && day % 10 == 0) {
// 			long noBytes = MemoryUtil.deepMemoryUsageOf(this);
// 			System.out.println("Total: " + noBytes);
// 			HostPopulation hp = demes.get(1);
// 			noBytes = MemoryUtil.deepMemoryUsageOf(hp);
// 			System.out.println("One host population: " + noBytes);
// 			Host h = hp.getRandomHostS();
// 			noBytes = MemoryUtil.deepMemoryUsageOf(h);
// 			System.out.println("One susceptible host with " +  h.getHistoryLength() + " previous infection: " + noBytes);
// 			//h.printHistory();
// 			if (getI() > 0) {
// 				Virus v = getRandomInfection();
// 				noBytes = MemoryUtil.memoryUsageOf(v);
// 				System.out.println("One virus: " + noBytes);
// 				noBytes = MemoryUtil.deepMemoryUsageOf(VirusTree.getTips());
// 				System.out.println("Virus tree: " + noBytes);
// 			}
// 		}
		
	}

	public void printHeader(PrintStream stream) {
		stream.print("date\tdiversity\ttmrca\tnetau\tserialInterval\tantigenicDiversity\ttotalN\ttotalS\ttotalI\ttotalT\ttotalCases\ttotalUnexposed");
		for (int i = 0; i < params.demeCount; i++) {
			HostPopulation hp = demes.get(i);
			hp.printHeader(stream);
		}
		stream.println();
	}
	
	public void printState(PrintStream stream) {
		stream.printf("%.4f\t%.4f\t%.4f\t%.4f\t%.5f\t%.4f\t%d\t%d\t%d\t%d\t%d\t%d", getDate(), getDiversity(), getTmrca(), getNetau(), getSerialInterval(), getAntigenicDiversity(), getN(), getS(), getI(), getT(), getCases(), getUnexposed());
		for (int i = 0; i < params.demeCount; i++) {
			HostPopulation hp = demes.get(i);
			hp.printState(stream);
		}
		stream.println();
	}	
	
	public void printVaccine() {
		Phenotype vaccineStrain = getLastVaccineStrain();
// 		try {
// 			sampleDb.beginTransaction(SqlJetTransactionMode.WRITE);
// 			ISqlJetTable table = sampleDb.getTable("vaccines");
// 			table.insert(nVaccines() - 1, lastVaccineDate,
// 				vaccineStrain.getTraitA(),
// 				vaccineStrain.getTraitB(),
// 				idealVaccineStrain.getTraitA(),
// 				idealVaccineStrain.getTraitB()
// 			);
// 			sampleDb.commit();
// 		}
// 		catch(SqlJetException e) {
// 			throw new RuntimeException(e);
// 		}

		vaccineStream.printf("%d\t%.4f\t%.4f\t%.4f\t%.4f\n",
				(int) day, vaccineStrain.getTraitA(), vaccineStrain.getTraitB(),
				idealVaccineStrain.getTraitA(), idealVaccineStrain.getTraitB());
	}
	
	public void printNVaccinations() {
// 		try {
// 			sampleDb.beginTransaction(SqlJetTransactionMode.WRITE);
// 			ISqlJetTable table = sampleDb.getTable("nVaccinations");
// 			for(int vaccineId = 0; vaccineId <= deployedVaccine; vaccineId++) {
// 				for(int demeId =  0; demeId < demes.size(); demeId++) {
// 					table.insert(getDate(), vaccineId, demeId, "S", demes.get(demeId).getNVaccinatedS(vaccineId));
// 					table.insert(getDate(), vaccineId, demeId, "I", demes.get(demeId).getNVaccinatedI(vaccineId)); 
// 				}
// 			}
// 			sampleDb.commit();
// 		}
// 		catch(SqlJetException e) {
// 			throw new RuntimeException(e);
// 		}
	}
	
	public void printSummary() {
	
		try {
			File summaryFile = new File("out.summary");
			summaryFile.delete();
			summaryFile.createNewFile();
			PrintStream summaryStream = new PrintStream(summaryFile);
			summaryStream.printf("parameter\tfull\n");
			summaryStream.printf("endDate\t%.4f\n", getDate());
			summaryStream.printf("diversity\t%.4f\n", mean(diversityList));
			summaryStream.printf("tmrca\t%.4f\n", mean(tmrcaList));
			summaryStream.printf("netau\t%.4f\n", mean(netauList));		
			summaryStream.printf("serialInterval\t%.5f\n", mean(serialIntervalList));	
			summaryStream.printf("antigenicDiversity\t%.4f\n", mean(antigenicDiversityList));	
			summaryStream.printf("N\t%.4f\n", mean(nList));		
			summaryStream.printf("S\t%.4f\n", mean(sList));		
			summaryStream.printf("I\t%.4f\n", mean(iList));		
			summaryStream.printf("R\t%.4f\n", mean(tList));		
			summaryStream.printf("cases\t%.4f\n", mean(casesList));	
			summaryStream.printf("Unexposed\t%.4f\n", mean(unexposedList));
			
			summaryStream.close();
		} catch(IOException ex) {
			System.out.println("Could not write to file"); 
			System.exit(0);
		}
	
	}
	
	private double mean(List<Double> list) {
		double mean = 0;
		if(!list.isEmpty()) {
			for (Double item : list) {
				mean += (double) item;
    		}
    	mean /= (double) list.size();
  		}
  		return mean;
	}	
		
	public void updateDiversity() {

		diversity = 0.0;
		tmrca = 0.0;
		antigenicDiversity = 0.0;		
		netau = 0.0;
		serialInterval = 0.0;
		
		double coalCount = 0.0;	
		double coalOpp = 0.0;
		double coalWindow = params.netauWindow / 365.0;
		int sampleCount = params.diversitySamplingCount;
		
		for (int i = 0; i < sampleCount; i++) {
			Virus vA = getRandomInfection();
			Virus vB = getRandomInfection();
			if (vA != null && vB != null) {
				double dist = vA.distance(vB);
				diversity += dist;
				if (dist > tmrca) {
					tmrca = dist;
				}
				antigenicDiversity += vA.antigenicDistance(vB);
				coalOpp += coalWindow;
				coalCount += vA.coalescence(vB, coalWindow);
				serialInterval += vA.serialInterval();
			}
		}	
	
		diversity /= (double) sampleCount;
		tmrca /= 2.0;	
		if ( tmrca > params.tmrcaLimit) {
		    System.out.println("Exceeding diversity threshold (TMRCA)." + tmrca);
		    String filename = "out.tmrcaLimit";
		    File excessFile = new File(filename);	
			try {
				excessFile.createNewFile();		    
			}
			catch (IOException ex) {
				System.out.println("Could not write to file"); 
				System.exit(0);
			}
		    System.exit(0);		
		}	
		netau = coalOpp / coalCount;
		serialInterval /= (double) sampleCount;
		antigenicDiversity /= (double) sampleCount;			
		
	}	
	
	public void pushLists() {
		diversityList.add(diversity);
		tmrcaList.add(tmrca);
		netauList.add(netau);
		serialIntervalList.add(serialInterval);
		antigenicDiversityList.add(antigenicDiversity);
		nList.add((double) getN());
		sList.add((double) getS());
		iList.add((double) getI());
		tList.add((double) getT());
		casesList.add((double) getCases());		
		unexposedList.add((double) getUnexposed());
	}
				
	public void resetCases() {
		for (int i = 0; i < params.demeCount; i++) {	
			HostPopulation hp = demes.get(i);
			hp.resetCases();
		}
	}
		
	public void stepForward() {
		// Update vaccine every vaccineStep
		// This happens even if vaccination is off, for the sake of comparison
		
		if(timestep % updateVaccineEvery == 0) {
			System.err.printf("Updating vaccine on day %s\n",day);
			updateVaccineStrain();
			System.err.printf("Most recent vaccine: %s\n",nVaccines()-1);
			System.err.printf("Currently deployed vaccine : %s\n",deployedVaccine);
			System.err.printf("Deploy new vaccine on day %s\n",lastVaccineDate*365+params.vaccineLag);
		}
		
		if(day == lastVaccineDate * 365.0 + params.vaccineLag){
			System.err.printf("Deploy new vaccine on day %s\n",day);
			System.err.printf("Previously deployed vaccine : %s\n",deployedVaccine);
			updateDeployedVaccine();
			System.err.printf("Updated deployed vaccine : %s\n",deployedVaccine);
		}
		for (int i = 0; i < params.demeCount; i++) {		
			HostPopulation hp = demes.get(i);
			
			// Only vaccinate individuals if the time of year falls within the vaccine window
			if(params.vaccinate) {
				if(day % 365.0 == params.deployDay){
					vaccineSeason = true;
					vaccineSeasonEndDay = day + params.vaccineWindow;
					System.err.printf("start vaccine season\n");
					System.err.printf("vaccine season ends on day %s\n",vaccineSeasonEndDay);
				}
				if(vaccineSeason){
					hp.vaccinate(deployedVaccine);
				}
			}
			hp.stepForward(day);
			for (int j = 0; j < params.demeCount; j++) {
				if (i != j) {
					HostPopulation hpOther = demes.get(j);
					hp.betweenDemeContact(hpOther, day);
				}
			}
		}
		
		timestep += 1;					
		day = timestep*params.deltaT;
		if(day == vaccineSeasonEndDay){
			vaccineSeason = false;
			System.err.printf("vaccine season ended on day %s\n",day);
		}
	}
	
	public void run() {
	
		try {
			
			// Summary time series
			File seriesFile = new File("out.timeseries");		
			seriesFile.delete();
			seriesFile.createNewFile();
			PrintStream seriesStream = new PrintStream(seriesFile);
			System.out.println("day\tdiversity\ttmrca\tnetau\tserialInterval\tantigenicDiversity\tN\tS\tI\tT\tcases\tUnexposed");
			printHeader(seriesStream);
			
			// Vaccine time series
			File vaccineFile = new File("vaccine.timeseries");
			vaccineFile.delete();
			vaccineFile.createNewFile();
			vaccineStream = new PrintStream(vaccineFile);
			vaccineStream.println("day\tactualA\tactualB\tidealA\tidealB");
			
			while (timestep <= endstep) {
									
				if (day % (double) params.printStep < params.deltaT) {			
					updateDiversity();
					printState();
					if (day > params.burnin) {
						printState(seriesStream);
						if(params.vaccinate && day > params.deployDay) { // only print vaccines after first vaccinations happen
							printNVaccinations();
						}
						pushLists();
					}
					resetCases();
				}
				
				if (getI()==0) {
					String exfilename = "out.extinct";
//					updateDiversity();
					printState(seriesStream);
//					updateVaccineStrain();
//					printVaccine();
					File extinctFile = new File(exfilename);
					try {
						extinctFile.createNewFile();
					}
					catch (IOException ex) {
						System.out.println("Could not write to file");
						System.exit(0);
					}
					System.exit(0);
					break;
				}
				
				stepForward();
			}
			seriesStream.close();
		} catch(IOException ex) {
			System.out.println("Could not write to file"); 
			System.exit(0);
		}	
		
		// tree reduction
		VirusTree.pruneTips();
		VirusTree.markTips();		
		VirusTree.reroot();		
	
		// tree prep
		makeTrunk();
		VirusTree.fillBackward();			
		VirusTree.sortChildrenByDescendants();
		VirusTree.setLayoutByDescendants();
		VirusTree.streamline();			
		
		// rotation
		if (params.pcaSamples) {
			VirusTree.rotate();
			VirusTree.flip();
		}
		
		
		// Summary
		printSummary();		
		VirusTree.printMKSummary();		// appends to out.summary
		
		//printHostPopulation();


		if (!params.reducedOutput) {	
		
			// tip and tree output	
			VirusTree.printTips();				
			VirusTree.printBranches();	
			VirusTree.printNewick();
					
			// immunity output
			VirusTree.updateRange();
			VirusTree.printRange();
			if (params.immunityReconstruction) {
				printImmunity();
			}
		
			// detailed output
//			if (params.detailedOutput) {
//				printHostPopulation();
//			}					
			
		}
		
	}
	
	public void initSampleDb() {
// 		try {
// 			File outputDbFile = new File("samples.sqlite");
// 			outputDbFile.delete();
// 			sampleDb = new SqlJetDb(outputDbFile, true);
// 			sampleDb.open();
// 			
// 			sampleDb.beginTransaction(SqlJetTransactionMode.WRITE);
// 			//sampleDb.createTable("CREATE TABLE hosts (hostId INTEGER);");
// 			sampleDb.createTable("CREATE TABLE vaccines (vaccineId INTEGER, date REAL, actualA REAL, actualB REAL, idealA REAL, idealB REAL);");
// 			//sampleDb.createTable("CREATE TABLE infections (hostId INTEGER, traitA REAL, traitB REAL);");
// 			//sampleDb.createTable("CREATE TABLE immuneHistory (hostId INTEGER, traitA REAL, traitB REAL);");
// 			//sampleDb.createTable("CREATE TABLE vaccineHistory (hostId INTEGER, vaccineId INTEGER);");
// 			sampleDb.createTable("CREATE TABLE nVaccinations (date REAL, vaccineId INTEGER, demeId INTEGER, state TEXT, nVaccinated INTEGER)");
// 			sampleDb.commit();
// 		}
// 		catch(SqlJetException e) {
// 			throw new RuntimeException(e);
// 		}
	}
	
	public double getDay() {
		return day;
	}
	
	// measured in years, starting at burnin
	public Double getDate() {
		return ((double) day - (double) params.burnin ) / 365.0;
	}
	
	public boolean dayIsInteger() {
   		return Math.ceil(day) - Math.floor(day) == 0;
	}		
}
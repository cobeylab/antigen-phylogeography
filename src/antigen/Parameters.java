package antigen;
/* Stores parameters for use across simulation */
/* Start with parameters in source, implement input file later */

public class Parameters {
	public Integer randomSeed; 	
	
	// simulation parameters
	public Integer burnin;
	public Integer endDay; 
	public Double deltaT;                                 	// number of days to move forward in a single timestep	
	public Integer printStep;									// print to out.timeseries every 10 days
	public Double tipSamplingRate;						// in samples per deme per day
	public Integer tipSamplesPerDeme;
	public Boolean tipSamplingProportional;				// whether to sample proportional to prevalance
	public Double treeProportion;							// proportion of tips to use in tree reconstruction
	public Integer diversitySamplingCount;					// how many samples to draw to calculate diversity, Ne*tau, serial interval
	public Integer netauWindow;								// window in days to calculate Ne*tau		
	public Boolean immunityReconstruction;				// whether to print immunity reconstruction to out.immunity
// 	public Boolean memoryProfiling;						// requires -javaagent:classmexer.jar to run
	public Double yearsFromMK;
	public Boolean pcaSamples;							// whether to rotate and flip virus tree
	public Boolean reducedOutput;						// whether to output only out.summary and out.timeseries
	public Double tmrcaLimit;
	
	public Boolean[] startAtEquilibriumInfected;
	public Boolean[] startAtEquilibriumImmune;
	
	// metapopulation parameters
	public Integer demeCount;
	public String[] demeNames;
	public int[] initialNs;	
	
	// host parameters
	public double[] birthRate;				// in births per individual per day, 1/30 years = 0.000091
	public double[]	deathRate;				// in deaths per individual per day, 1/30 years = 0.000091
	public Boolean swapDemography;				// whether to keep overall population size constant
		
	// epidemiological parameters
	public int[] initialIs;							// in individuals
	public int[] initialRs; // ONLY USED FOR OUTPUT PARAMETERS
	public Integer initialDeme;						// index of deme where infection starts, 1..n
	public Double initialPrR; 					// as proportion of population
	public Double beta; // 0.3					// in contacts per individual per day
	public Double nu; //0.2						// in recoveries per individual per day
	public Double betweenDemePro;				// relative to within-deme beta	
	public Double[][] contactMatrix;            // relative to betweenDemePro
	
	// transcendental immunity
	public Boolean transcendental;
	public Double immunityLoss;					// in R->S per individual per day
	public Double initialPrT;

	// seasonal betas
	public double[] demeBaselines;
	public double[] demeAmplitudes;
	public double[] demeOffsets;				// relative to the year
	
	// phenotype parameters
	public Double muPhenotype; 					// in mutations per individual per day

	// parameters specific to GeometricPhenotype
	public Double smithConversion;					// multiplier to distance to give cross-immunity
	public Double homologousImmunity;				// immunity raised to antigenically identical virus
	public Double initialTraitA;	
	public Double meanStep; 
	public Double sdStep; 
	public Boolean mut2D;						// whether to mutate in a full 360 degree arc
	
	// vaccination parameters
	public Boolean vaccinate;					// whether to include vaccination
	public Double vaccineStep;					// how many days between re-calculation of ideal vaccine
	public Double vaccineSD;						// SD of error relative to ideal vaccine
	public Double vaccineImmuneBreadth;			// relative immune protection provided by vaccine relative to infection:
																// 1.0 -> same protection; 2.0 -> double protection (risk halved)
	public Double[] vaccinationRate; 							// Average time between vaccinations for 0.003 is a little less than once a year
																// : (1/0.003) = 333.3 days
	public Double vaccineLag;					// number of days between vaccine selection and deployment (needs to be less than 1 year) 
	public Boolean synchronizeVaccination; // whether to synchronize vaccination for 4 months
	public Double deployDay; //how many days between start of vaccine 'season'
	public Double vaccineWindow; //how many days the vaccine 'season' lasts
	
	// Periodic sampling of hosts into SQLite db
	// TODO: not implemented
	public Double sampleStep;
	public Integer sampleHostsPerDeme;
}
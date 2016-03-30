package antigen;
/* Holds random number genator necessities */
/* Trying to encapsulate this, so the RNG particulars can be changed if necessary */ 
/* Completely static class, allows no instances to be instantiated */

//import cern.jet.random.*;
import cern.jet.random.*;
import cern.jet.random.engine.*;

public class Random {
	static RandomEngine rng;
	static Uniform uniform;
	static Normal normal;
	static Gamma gamma;
	static Poisson poisson;
	
	public static int seed() {
		RandomEngine tmpRng = RandomEngine.makeDefault();
		int seed = tmpRng.nextInt();
		rng = new MersenneTwister(seed);
		initDistributions();
		return seed;
	}
	
	public static void seed(int seed) {
		rng = new MersenneTwister(seed);
		initDistributions();
	}
	
	private static void initDistributions() {
		uniform = new Uniform(rng);
		normal = new Normal(0.0, 1.0, rng);
		gamma = new Gamma(1.0, 1.0, rng);
		poisson = new Poisson(1.0, rng);
	}
	
	public static int nextInt(int from, int to) {
		return uniform.nextIntFromTo(from, to);
	}	
	
	public static double nextDouble() {
		return uniform.nextDouble();		
	}
	
	public static double nextDouble(double from, double to) {
		return uniform.nextDoubleFromTo(from, to);		
	}	

	public static double nextNormal() {
		return normal.nextDouble();
	}
	
	public static double nextNormal(double mean, double sd) {
		return normal.nextDouble(mean, sd);
	}
	
	// tuned with alpha and beta, matching Mathematica's notation
	public static double nextGamma(double alpha, double beta) {
		return gamma.nextDouble(alpha, 1/beta);
	}	
	
	public static int nextPoisson(double lambda) {
		return poisson.nextInt(lambda);
	}
	
	public static boolean nextBoolean(double p) {
		boolean x = false;
		if (nextDouble() < p) {
			x = true;
		}
		return x;
	}
}

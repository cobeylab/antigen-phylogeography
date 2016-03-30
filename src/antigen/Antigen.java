package antigen;
import java.io.FileReader;

import com.google.gson.Gson;

/* Implements an individual-based model in which the infection's genealogical history is tracked through time */

class Antigen {
    public static void main(String[] args) {
    	// Load parameters from parameters.json
		Parameters params;
		try {
			String filename;
			if(args.length > 0) {
				filename = args[0];
			}
			else {
				filename = "parameters.json";
			}
			params = new Gson().fromJson(new FileReader(filename), Parameters.class);
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
		
		// Run simulation
		Simulation sim = new Simulation(params);
		sim.run();	
	}
}


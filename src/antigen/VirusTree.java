package antigen;
/* Stores a list of Viruses that have sampled during the course of the simulation */

import java.util.*;
import java.io.*;

public class VirusTree {

	// fields
	public static Simulation sim;
	public static Parameters params;
	public static Virus root;	
	private static List<Virus> tips = new ArrayList<Virus>();
	
	public static double xMin;
	public static double xMax;
	public static double yMin;
	public static double yMax;
	public static double zMin;
	public static double zMax;	
	
	static final Comparator<Virus> descendantOrder = new Comparator<Virus>() {
		public int compare(Virus v1, Virus v2) {
			Integer descendantsV1 = new Integer(getNumberOfDescendants(v1));
			Integer descendantsV2 = new Integer(getNumberOfDescendants(v2));
			return descendantsV1.compareTo(descendantsV2);
		}
	};	
		
	// static methods
	public static void add(Virus v) {		
		tips.add(v);
	}
	public static void clear() {
		tips.clear();
	}
	public static List<Virus> getTips() {
		return tips;
	}
	public static Virus getRoot() {
		return root;
	}
	
	// go through tips and find TMRCA
	public static Virus getTMRCA() {
		if(tips.size() == 0) {
			return null;
		}
		
		Virus tmrca = tips.get(0);
		for(int i = 1; i < tips.size(); i++) {	
			tmrca = tmrca.commonAncestor(tips.get(i));
		}
		return tmrca;
	}
	
	// reroot tree at TMRCA rather than urVirus
	public static void reroot() {
		root = getTMRCA();
	}
	
	// return a random tip that lies between year from and year to
	public static Virus getRandomTipFromTo(double from, double to) {
	
		// fill temporary list
		List<Virus> select = new ArrayList<Virus>();
		for (Virus v : tips) {
			double x = v.getBirth();
			if (x >= from && x < to) {
				select.add(v);
			}
		}
		
		// pull random virus from this list
		Virus rV = null;
		if (select.size() > 0) {	
			int index = Random.nextInt(0,select.size()-1);
			rV = select.get(index);
		}
		return rV;
		
	}
	
	public static int getDemeCount(int d) {
		int count = 0;
		for (Virus v : tips) {
			if (v.getDeme() == d) {
				count++;
			}
		}
		return count;
	}	
		
	// work backwards for each sample filling the children lists
	public static void fillBackward() {
	
		for (Virus child : tips) {
			Virus parent = child.getParent();
			while (parent != null) {
				parent.addChild(child);
				parent.incrementCoverage();
				child = parent;
				parent = child.getParent();
			}
		}
	
	}
	
	public static void dropTips() {
	
		List<Virus> reducedTips = new ArrayList<Virus>();
		for (Virus v : tips) {
			if (Random.nextBoolean(params.treeProportion)) {
				reducedTips.add(v);
			}
		}
		tips = reducedTips;
	
	}

	// marking to by time, not proportional to prevalence
	public static void markTips() {
	
//		for (Virus v : tips) {
//			if (Random.nextBoolean(params.treeProportion)) {
//				while (v.getParent() != null) {
//					v.mark();
//					v = v.getParent();
//				}
//			}
//		}
		
		for (double i = 0; i < sim.getDate(); i+=0.1) {
			Virus v = getRandomTipFromTo(i,i+0.1);
			if (v != null) {
				while (v.getParent() != null) {
					v.mark();
					v = v.getParent();
				}
			}
		}
		
	}
		
	// prune tips
	public static void pruneTips() {
	
		List<Virus> reducedTips = new ArrayList<Virus>();
		for (int d = 0; d < params.demeCount; d++) {
			double keepProportion = (double) params.tipSamplesPerDeme / (double) getDemeCount(d);
			for (Virus v : tips) {
				if (Random.nextBoolean(keepProportion) && v.getDeme() == d) {
					reducedTips.add(v);
				}
			}
		}
		tips = reducedTips;
	
	}
	
	// returns virus v and all its descendents via a depth-first traversal
	public static List<Virus> postOrderNodes(Virus v) {
		List<Virus> vNodes = new ArrayList<Virus>();
		vNodes.add(v);
		vNodes = postOrderChildren(vNodes);
		return vNodes;
	}
	
	public static List<Virus> postOrderNodes() {
		return postOrderNodes(root);
	}	
	
	// returns virus v and all its descendents via a depth-first traversal
	public static List<Virus> postOrderChildren(List<Virus> vNodes) {
	
		Virus last = vNodes.get(vNodes.size()-1);
	
		for (Virus child : last.getChildren()) {
			vNodes.add(child);
			postOrderChildren(vNodes);
		}
		
		return vNodes;
	
	}


	// Count total descendents of a Virus, working through its children and its children's children
	public static int getNumberOfDescendants(Virus v) {
	
		int numberOfDescendants = v.getNumberOfChildren();
		
		for (Virus child : v.getChildren()) {
			numberOfDescendants += getNumberOfDescendants(child);
		}
		
		return numberOfDescendants;
		
	}
	
	public static int getNumberOfDescendants() {
		return getNumberOfDescendants(root);
	}
		
	// sorts children lists so that first member is child with more descendents than second member
	public static void sortChildrenByDescendants(Virus v) {
		
		List<Virus> children = v.getChildren();
		Collections.sort(children, descendantOrder);
		
		for (Virus child : children) {
			sortChildrenByDescendants(child);
		}
				
	}	
	
	public static void sortChildrenByDescendants() {
		sortChildrenByDescendants(root);
	}
	
	// sets Virus layout based on a postorder traversal
	public static void setLayoutByDescendants() {
	
		List<Virus> vNodes = postOrderNodes();
		
		// set layout of tips based on traversal
		double y = 0;
		for (Virus v : vNodes) {
//			if (tips.contains(v)) {
			if (v.isTip()) {
				v.setLayout(y);
				y++;
			}
		}
		
		// update layout of internal nodes
		Collections.reverse(vNodes);
		for (Virus v : vNodes) {
			if (v.getNumberOfChildren() > 0) {
				double mean = 0;
				for (Virus child : v.getChildren()) {
					mean += child.getLayout();
				}
				mean /= v.getNumberOfChildren();
				v.setLayout(mean);
			}
		}
		
	}	
	
	// looks at a virus and its grandparent, if traits are identical and there is no branching
	// then make virus child rather than grandchild
	// returns v.parent after all is said and done
	public static Virus collapse(Virus v) {
	
		Virus vp = null;
		Virus vgp = null;
		if (v.getParent() != null) {
			vp = v.getParent();
			if (vp.getParent() != null) {
				vgp = vp.getParent();
			}
		}

		if (vp != null && vgp != null) {
	//		if (vp.getNumberOfChildren() == 1 && v.getPhenotype() == vp.getPhenotype() && v.isTrunk() == vp.isTrunk() && v.getDeme() == vp.getDeme()) {
		
			if (vp.getNumberOfChildren() == 1) {
		
				List<Virus> vgpChildren = vgp.getChildren();
				int vpIndex =  vgpChildren.indexOf(vp);
				
				if (vpIndex >= 0) {
				
					// replace virus as child of grandparent
					vgpChildren.set(vpIndex, v);
				
					// replace grandparent as parent of virus
					v.setParent(vgp);
				
					// erase parent
					vp = null;
				
				}
		
			}
		}
		
		return v.getParent();

	}
	
	// walks backward using the list of tips, collapsing where possible
	public static void streamline() {
		
		for (Virus v : tips) {
			Virus vp = v;
			while (vp != null) {
				vp = collapse(vp);
			}
		}
		
	}
	
	// rotate the 2d euclidean space using PCA, returning an x-axis with maximum variance
	public static void rotate() {
			
		// load a 2d array with phenotypes
		
		List<Virus> virusList = postOrderNodes();
		int n = virusList.size();
		int m = 2;
		
		double[][] input = new double[n][m];
		
		for (int i = 0; i < n; i++) {
			Virus v = virusList.get(i);
			Phenotype p = v.getPhenotype();
			double x = p.getTraitA();
			double y = p.getTraitB();	
			input[i][0] = x;
			input[i][1] = y;				
		}
		
		// project this array
		
		double[][] projected = SimplePCA.project(input);
		
		// reset phenotypes based on projection
		
		for (int i = 0; i < n; i++) {
			Virus v = virusList.get(i);
			Phenotype p = v.getPhenotype();
			double x = projected[i][0];
			double y = projected[i][1];				
			p.setTraitA(x);
			p.setTraitB(y);					
		}
	}
	
	// flips the 2d euclidean space so that first sample is always to the left of the last sample
	public static void flip() {

		List<Virus> virusList = postOrderNodes();
		int n = virusList.size();	
		
		// find first and last virus			
		Virus firstVirus = virusList.get(0);
		Virus lastVirus = virusList.get(0);
		double firstDate = firstVirus.getBirth();
		double lastDate = lastVirus.getBirth();
				
		for (Virus v : virusList) {
			if (v.getBirth() < firstDate) {
				firstDate = v.getBirth();
				firstVirus = v;
			}
			if (v.getBirth() > lastDate) {
				lastDate = v.getBirth();
				lastVirus = v;
			}				
		}
		
		// is the x-value of first virus greater than the x-value of last virus?
		// if so, flip
		
		Phenotype p = firstVirus.getPhenotype();
		double firstX = p.getTraitA();
		p = lastVirus.getPhenotype();
		double lastX = p.getTraitA();		
		
		if (firstX > lastX) {
		
			// I think that postOrderNodes() has replicates in it, need to go through some hoops because of this
			double[] input = new double[n];
		
			for (int i = 0; i < n; i++) {
				Virus v = virusList.get(i);
				p = v.getPhenotype();		
				input[i] = p.getTraitA();;
			}
			
			for (int i = 0; i < n; i++) {
				Virus v = virusList.get(i);
				p = v.getPhenotype();
				double x = -1*input[i];		
				p.setTraitA(x);
			}				
		
		}
	
	}
	
	// walks through list of nodes and update min and max ranges appropriately
	public static void updateRange() {
	
		xMin = 0.0;
		xMax = 0.0;
		yMin = 0.0;
		yMax = 0.0;
		zMin = 0.0;
		zMax = 0.0;		
	
		for (Virus v : postOrderNodes()) {
		
			Phenotype p =  v.getPhenotype();
			double x = p.getTraitA();
			double y = p.getTraitB();
			if (xMin > x) { xMin = x; }
			if (xMax < x) { xMax = x; }
			if (yMin > y) { yMin = y; }
			if (yMax < y) { yMax = y; }	
		
		}		
		
		xMin = Math.floor(xMin) - 10;
		xMax = Math.ceil(xMax) + 10;
		yMin = Math.floor(yMin) - 10;
		yMax = Math.ceil(yMax) + 10;
		zMin = Math.floor(zMin) - 10;
		zMax = Math.ceil(zMax) + 10;		
	
	}

	public static void printRange() {
		
		try {
			File rangeFile = new File("out.range");
			rangeFile.delete();
			rangeFile.createNewFile();
			PrintStream rangeStream = new PrintStream(rangeFile);
			rangeStream.printf("%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n", xMin, xMax, yMin, yMax, zMin, zMax);
			rangeStream.close();
		} catch(IOException ex) {
			System.out.println("Could not write to file"); 
			System.exit(0);
		}
		
	}
	
	public static void printTips() {
		
		try {
			File tipFile = new File("out.tips");
			tipFile.delete();
			tipFile.createNewFile();
			PrintStream tipStream = new PrintStream(tipFile);
			tipStream.printf("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\"\n", "name", "year", "trunk", "tip", "mark", "location", "layout", "ag1", "ag2");
			for (int i = 0; i < tips.size(); i++) {
				Virus v = tips.get(i);			
				tipStream.printf("\"%s\",%.4f,%d,%d,%d,%d,%.4f,%s\n", v, v.getBirth(), v.isTrunk()?1:0, v.isTip()?1:0, v.isMarked()?1:0, v.getDeme(), v.getLayout(), v.getPhenotype());
			}
			tipStream.close();
		} catch(IOException ex) {
			System.out.println("Could not write to file"); 
			System.exit(0);
		}
		
	}
	
	public static void printBranches() {
		
		try {
			File branchFile = new File("out.branches");
			branchFile.delete();
			branchFile.createNewFile();
			PrintStream branchStream = new PrintStream(branchFile);
			for (Virus v : postOrderNodes()) {
				if (v.getParent() != null) {
					Virus vp = v.getParent();
					branchStream.printf("{\"%s\",%.4f,%d,%d,%d,%d,%.4f,%s}\t", v, v.getBirth(), v.isTrunk()?1:0, v.isTip()?1:0, v.isMarked()?1:0, v.getDeme(), v.getLayout(), v.getPhenotype());
					branchStream.printf("{\"%s\",%.4f,%d,%d,%d,%d,%.4f,%s}\t", vp, vp.getBirth(), vp.isTrunk()?1:0, vp.isTip()?1:0, v.isMarked()?1:0, vp.getDeme(), vp.getLayout(), vp.getPhenotype());
					branchStream.printf("%d\n", vp.getCoverage());
				}
			}
			branchStream.close();
		} catch(IOException ex) {
			System.out.println("Could not write to file"); 
			System.exit(0);
		}
		
	}
	
	// assess node in building Newick string
	public static Virus assessNode(Virus v, List<Virus> visited, PrintStream treeStream) {
	
		Virus returnVirus = null;
		boolean printHeight = false;
	
		// if virus has multiple children, return first child that has not been visited
		if (v.getNumberOfChildren() > 1) {
			boolean childrenVisited = true;
			for (int i = 0; i < v.getNumberOfChildren(); i++) {
				Virus vc = v.getChildren().get(i);
				if (!visited.contains(vc)) {
					if (i == 0) {
						treeStream.print("(");
					}
					else {
						treeStream.print(",");
					}
					childrenVisited = false;
					returnVirus = vc;
					break;
				}
			}
			// failure, all children visited, return to parent
			if (childrenVisited) {
				treeStream.print(")");	
				printHeight = true;
				returnVirus = v.getParent();
			}
		}
		
		// if tip is encountered, print tip, return to parent
		if (v.getNumberOfChildren() == 0) {
			treeStream.print(v.toString());
			printHeight = true;
			returnVirus = v.getParent();		
		}			
		
		// walk down (or up) branches 
		if (v.getNumberOfChildren() == 1) {
			Virus vc = v.getChildren().get(0);
			if (!visited.contains(vc)) {
				returnVirus = vc;
			}
			else {
				returnVirus = v.getParent();
			}
		}
		
		// find height, walk back until a parent with a split occurs
		if (printHeight && v.getParent() != null) {

			treeStream.printf("[&antigenic={%s}]", v.getPhenotype() );
		
			Virus vp = v.getParent();
			while (vp.getNumberOfChildren() == 1 && vp.getParent() != null) {
				vp = vp.getParent();
			}
			double height = v.getBirth() - vp.getBirth();
			treeStream.printf(":%.4f", height);	

		}
		
		return returnVirus;
	
	}
	
	public static void printNewick() {
	
		try {
			File treeFile = new File("out.trees");
			treeFile.delete();
			treeFile.createNewFile();
			PrintStream treeStream = new PrintStream(treeFile);
			
			List<Virus> visited = new ArrayList<Virus>();
				
			// start at root
			Virus v = root;
			visited.add(v);			
			
			while (v != null) {
				
				v = assessNode(v, visited, treeStream);		
				visited.add(v);
							
			}
			treeStream.println();
			
			treeStream.close();
		} catch(IOException ex) {
			System.out.println("Could not write to file"); 
			System.exit(0);
		}
	
	}
	
	public static int sideBranchMutations() {
		int count = 0;
		for (Virus v : postOrderNodes()) {
			if (v.getParent() != null && v.getBirth() < sim.getDate() - params.yearsFromMK) {
				Virus vp = v.getParent();
				if (!v.isTrunk() && !vp.isTrunk() && v.getPhenotype() != vp.getPhenotype()) {
					count++;
				}
			}
		}
		return count;
	}	
	
	public static double sideBranchOpportunity() {
		double time = 0;
		for (Virus v : postOrderNodes()) {
			if (v.getParent() != null && v.getBirth() < sim.getDate() - params.yearsFromMK) {
				Virus vp = v.getParent();
				if (!v.isTrunk() && !vp.isTrunk()) {
					time += v.getBirth() - vp.getBirth();
				}
			}
		}
		return time;
	}	
	
	public static int trunkMutations() {
		int count = 0;
		for (Virus v : postOrderNodes()) {
			if (v.getParent() != null && v.getBirth() < sim.getDate() - params.yearsFromMK) {
				Virus vp = v.getParent();
				if (v.isTrunk() && vp.isTrunk() && v.getPhenotype() != vp.getPhenotype()) {
					count++;
				}
			}
		}
		return count;
	}	
	
	public static double trunkOpportunity() {
		double time = 0;
		for (Virus v : postOrderNodes()) {
			if (v.getParent() != null && v.getBirth() < sim.getDate() - params.yearsFromMK) {
				Virus vp = v.getParent();
				if (v.isTrunk() && vp.isTrunk()) {
					time += v.getBirth() - vp.getBirth();
				}
			}
		}
		return time;
	}		
	
	public static void printMKSummary() {
		
		try {
			PrintStream summaryStream = new PrintStream(new FileOutputStream("out.summary", true)); // append
			double sideBranchMut = (double) sideBranchMutations();
			double sideBranchOpp = sideBranchOpportunity();
			double sideBranchRate = sideBranchMut / sideBranchOpp;
			double trunkMut = (double) trunkMutations();
			double trunkOpp = trunkOpportunity();	
			double trunkRate = trunkMut / trunkOpp;		
			double mkRatio = trunkRate / sideBranchRate;
			summaryStream.printf("sideBranchRate\t%.4f\n", sideBranchRate);	
			summaryStream.printf("trunkRate\t%.4f\n", trunkRate);	
			summaryStream.printf("mkRatio\t%.4f\n", mkRatio);	
			summaryStream.close();
		} catch(IOException ex) {
			System.out.println("Could not write to file"); 
			System.exit(0);
		}
		
	}	
		
}
package antigen;

import java.util.*;

@SuppressWarnings("serial")
public class IterableBitSet extends BitSet implements Iterable<Integer> {
	
	public IterableBitSet() {
	}
	
	public IterableBitSet(int n) {
		super(n);
	}
	
	public int[] getSetBits() {
		int[] setBits = new int[cardinality()];
		int j = 0;
		for(int i = nextSetBit(0); i >= 0; i = nextSetBit(i + 1)) {
			setBits[j++] = i;
		}
		return setBits;
	}
	
	@Override
	public Iterator<Integer> iterator() {
		return new BitSetIterator(this);
	}
	
	@Override
	public Object clone() {
		IterableBitSet bs = new IterableBitSet();
		for(int i : this) {
			bs.set(i);
		}
		return bs;
	}
}

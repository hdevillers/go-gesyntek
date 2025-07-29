package kmer

import (
	"testing"
)

// Test no split required
func TestNoSplitRequired(t *testing.T) {
	kd := NewKSplit(5)
	seq := []byte("ACGCTCGCGCGATCGATCGAGCTATGCGTC") // 30 bases

	err := kd.SplitSeq(&seq)

	if err != nil {
		t.Errorf("Unexpected error occurred while splitting a sequence: %s", err.Error())
	}

	if kd.NSkipped != 0 {
		t.Errorf("Sequence splitting process skipped %d bases as unsupported nucleotides while no bases should have been skipped.", kd.NSkipped)
	}

	if kd.NTooShort != 0 {
		t.Errorf("Sequence splitting process skipped %d bases in too short fragment(s) while no bases should have been skipped.", kd.NSkipped)
	}

	if len(kd.SeqSplit) != 1 {
		t.Errorf("Expected no splits but found %d split(s).", len(kd.SeqSplit)-1)
	}

	if len(kd.SeqSplit[0]) != len(seq) {
		t.Errorf("Expected kept sequence length of %d but found %d.", len(seq), len(kd.SeqSplit[0]))
	}

	s1 := string(seq)
	s2 := string(kd.SeqSplit[0])
	if s1 != s2 {
		t.Errorf("Splitting process has change nucleotide from %s to %s.", s1, s2)
	}
}

// Test mixture of upper and lower cases
func TestUpperLowerCases(t *testing.T) {
	kd := NewKSplit(5)
	seq := []byte("ACGCtCGaGCgtTCGATCggGCTATGaGtc") // 30 bases

	err := kd.SplitSeq(&seq)

	if err != nil {
		t.Errorf("Unexpected error occurred while splitting a sequence: %s", err.Error())
	}

	if kd.NSkipped != 0 {
		t.Errorf("Sequence splitting process skipped %d bases as unsupported nucleotides while no bases should have been skipped.", kd.NSkipped)
	}

	if kd.NTooShort != 0 {
		t.Errorf("Sequence splitting process skipped %d bases in too short fragment(s) while no bases should have been skipped.", kd.NSkipped)
	}

	if len(kd.SeqSplit) != 1 {
		t.Errorf("Expected no splits but found %d split(s).", len(kd.SeqSplit)-1)
	}

	if len(kd.SeqSplit[0]) != len(seq) {
		t.Errorf("Expected kept sequence length of %d but found %d.", len(seq), len(kd.SeqSplit[0]))
	}

	s1 := string(seq)
	s2 := string(kd.SeqSplit[0])
	if s1 != s2 {
		t.Errorf("Splitting process has change nucleotide from %s to %s.", s1, s2)
	}
}

// Test one split for one base
func TestOneSimpleSplit(t *testing.T) {
	kd := NewKSplit(5)
	seq := []byte("ACGCTCGCGCGATCGNTCGAGCTATGCGTC") // 30 bases

	err := kd.SplitSeq(&seq)

	if err != nil {
		t.Errorf("Unexpected error occurred while splitting a sequence: %s", err.Error())
	}

	if kd.NSkipped != 1 {
		t.Errorf("Sequence splitting process skipped %d bases as unsupported nucleotides while 1 base should have been skipped.", kd.NSkipped)
	}

	if kd.NTooShort != 0 {
		t.Errorf("Sequence splitting process skipped %d bases in too short fragment(s) while no bases should have been skipped.", kd.NSkipped)
	}

	if len(kd.SeqSplit) != 2 {
		t.Errorf("Expected one split but found %d split(s).", len(kd.SeqSplit)-1)
	}

	if len(kd.SeqSplit[0]) != 15 {
		t.Errorf("Expected kept sequence length for the first split of 15 but found %d.", len(kd.SeqSplit[0]))
	}

	if len(kd.SeqSplit[1]) != 14 {
		t.Errorf("Expected kept sequence length for the second split of 14 but found %d.", len(kd.SeqSplit[1]))
	}
}

// Test one split for one base
func TestTwoSplits(t *testing.T) {
	kd := NewKSplit(5)
	seq := []byte("ACGCTCGXGCGATCGNTCGAGCTATGCGTC") // 30 bases

	err := kd.SplitSeq(&seq)

	if err != nil {
		t.Errorf("Unexpected error occurred while splitting a sequence: %s", err.Error())
	}

	if kd.NSkipped != 2 {
		t.Errorf("Sequence splitting process skipped %d bases as unsupported nucleotides while 2 bases should have been skipped.", kd.NSkipped)
	}

	if kd.NTooShort != 0 {
		t.Errorf("Sequence splitting process skipped %d bases in too short fragment(s) while no bases should have been skipped.", kd.NSkipped)
	}

	if len(kd.SeqSplit) != 3 {
		t.Errorf("Expected one split but found %d split(s).", len(kd.SeqSplit)-1)
	}

	if len(kd.SeqSplit[0]) != 7 {
		t.Errorf("Expected kept sequence length for the first split of 7 but found %d.", len(kd.SeqSplit[0]))
	}

	if len(kd.SeqSplit[1]) != 7 {
		t.Errorf("Expected kept sequence length for the second split of 7 but found %d.", len(kd.SeqSplit[1]))
	}

	if len(kd.SeqSplit[2]) != 14 {
		t.Errorf("Expected kept sequence length for the second split of 14 but found %d.", len(kd.SeqSplit[2]))
	}
}

// Test skip a fragment that is too short at the end
func TestTooShortFragmentEnd(t *testing.T) {
	kd := NewKSplit(5)
	seq := []byte("ACGCTCGCGCGATCGATCGAGCTATNCGTC") // 30 bases

	err := kd.SplitSeq(&seq)

	if err != nil {
		t.Errorf("Unexpected error occurred while splitting a sequence: %s", err.Error())
	}

	if kd.NSkipped != 1 {
		t.Errorf("Sequence splitting process skipped %d bases as unsupported nucleotides while 1 base should have been skipped.", kd.NSkipped)
	}

	if kd.NTooShort != 4 {
		t.Errorf("Sequence splitting process skipped %d bases in too short fragment(s) while 4 bases should have been skipped.", kd.NSkipped)
	}

	if len(kd.SeqSplit) != 1 {
		t.Errorf("Expected no splits but found %d split(s).", len(kd.SeqSplit)-1)
	}

	if len(kd.SeqSplit[0]) != 25 {
		t.Errorf("Expected kept sequence length of 25 but found %d.", len(kd.SeqSplit[0]))
	}
}

// Test skip a fragment that is too short at the start
func TestTooShortFragmentStart(t *testing.T) {
	kd := NewKSplit(5)
	seq := []byte("ACGCNCGCGCGATCGATCGAGCTATTCGTC") // 30 bases

	err := kd.SplitSeq(&seq)

	if err != nil {
		t.Errorf("Unexpected error occurred while splitting a sequence: %s", err.Error())
	}

	if kd.NSkipped != 1 {
		t.Errorf("Sequence splitting process skipped %d bases as unsupported nucleotides while 1 base should have been skipped.", kd.NSkipped)
	}

	if kd.NTooShort != 4 {
		t.Errorf("Sequence splitting process skipped %d bases in too short fragment(s) while 4 bases should have been skipped.", kd.NSkipped)
	}

	if len(kd.SeqSplit) != 1 {
		t.Errorf("Expected no splits but found %d split(s).", len(kd.SeqSplit)-1)
	}

	if len(kd.SeqSplit[0]) != 25 {
		t.Errorf("Expected kept sequence length of 25 but found %d.", len(kd.SeqSplit[0]))
	}
}

// Test skip a fragment in the middle
func TestTooShortFragmentMiddle(t *testing.T) {
	kd := NewKSplit(5)
	seq := []byte("ACGCTCGCGCGNTCGNTCGAGCTATGCGTC") // 30 bases

	err := kd.SplitSeq(&seq)

	if err != nil {
		t.Errorf("Unexpected error occurred while splitting a sequence: %s", err.Error())
	}

	if kd.NSkipped != 2 {
		t.Errorf("Sequence splitting process skipped %d bases as unsupported nucleotides while 2 bases should have been skipped.", kd.NSkipped)
	}

	if kd.NTooShort != 3 {
		t.Errorf("Sequence splitting process skipped %d bases in too short fragment(s) while no bases should have been skipped.", kd.NSkipped)
	}

	if len(kd.SeqSplit) != 2 {
		t.Errorf("Expected one split but found %d split(s).", len(kd.SeqSplit)-1)
	}

	if len(kd.SeqSplit[0]) != 11 {
		t.Errorf("Expected kept sequence length for the first split of 15 but found %d.", len(kd.SeqSplit[0]))
	}

	if len(kd.SeqSplit[1]) != 14 {
		t.Errorf("Expected kept sequence length for the second split of 14 but found %d.", len(kd.SeqSplit[1]))
	}
}

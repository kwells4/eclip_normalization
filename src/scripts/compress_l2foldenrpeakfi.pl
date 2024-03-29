use warnings;
use strict;

my $hashing_value = 100000;
my $fi = $ARGV[0];
my $output_fi = $fi.".compressed.bed";
open(O,">$output_fi");

my %peaks2l10p;
my %peaks2start;
my %read_hash;
my %peak_hash;
&readfi($fi);

my %overlap_hash;
for my $chr (keys %read_hash) {
    for my $str ("+","-") {

	my %deleted_peaks;
	my %kept_peaks;

	my @sorted_peaks = sort {$peaks2l10p{$chr}{$str}{$fi}{$b} <=> $peaks2l10p{$chr}{$str}{$fi}{$a}} keys %{$peaks2l10p{$chr}{$str}{$fi}};
	
	my $i=0;

	for my $peak1 (@sorted_peaks) {
	    my $verbose_flag = 0;

## now take any peaks that overlap and merge them
	    next if (exists $deleted_peaks{$peak1});
	    my ($p1chr,$p1pos,$p1str,$p1vsinput_l10p,$p1vsinput_l2fenr) = split(/\:/,$peak1);
	    my ($p1start,$p1stop) = split(/\-/,$p1pos);

	    my $p1x = int($p1start / $hashing_value);
	    my $p1y = int( $p1stop / $hashing_value);

	    for my $p1i ($p1x..$p1y) {
		for my $tocomp_peak (@{$read_hash{$chr}{$str}{$fi}{$p1i}}) {
		    print STDERR "comparing $peak1 $tocomp_peak\n" if ($verbose_flag == 1);
		    next if (exists $deleted_peaks{$tocomp_peak});
		    next if ($tocomp_peak eq $peak1);

		    my ($p2compchr,$p2comppos,$p2compstr,$p2compvsinput_l10p,$p2compvsinput_l2fenr) = split(/\:/,$tocomp_peak);
		    my ($p2compstart,$p2compstop) = split(/\-/,$p2comppos);

		    next if ($p2compstop <= $p1start);
		    next if ($p1stop <= $p2compstart);

		    #peak2 overlaps with peak1 and has a lower l10pval - remove it!
		    if ($p1vsinput_l10p >= $p2compvsinput_l10p) {
			print STDERR "discarding $tocomp_peak vs $peak1\n" if ($verbose_flag == 1);
			$deleted_peaks{$tocomp_peak} = 1;
		    } elsif ($p1vsinput_l10p < $p2compvsinput_l10p) {
			$deleted_peaks{$peak1} = 1;
			print STDERR "discarding $peak1 vs $tocomp_peak\n" if ($verbose_flag == 1);
		    } else {
			print STDERR "weird error shouldn't happen $peak1\n";
		    }
		}
	    }
	}
	
	for my $peak (@sorted_peaks) {
	    next if (exists $deleted_peaks{$peak});

	    my ($p1chr,$p1pos,$p1str,$p1vsinput_l10p,$p1vsinput_l2fenr) = split(/\:/,$peak);
            my ($p1start,$p1stop) = split(/\-/,$p1pos);
	    print O "$p1chr\t$p1start\t$p1stop\t$p1vsinput_l10p\t$p1vsinput_l2fenr\t$p1str\n";

	}
    }
}
close(O);


sub min {
    my $x = shift;
    my $y = shift;
    
    if ($x < $y) {
	return($x);
    } else {
	return($y);
    }
}

sub max {
    my $x = shift;
    my $y = shift;

    if ($x > $y) {
	return($x);
    } else {
	return($y);
    }
}

sub readfi {
    my $fi = shift;
    open(F,$fi);
    for my $line (<F>) {
	chomp($line);

	my @tmp = split(/\t/,$line);

	my $chr = $tmp[0];
	my $str = $tmp[5];
	my $start = $tmp[1];
	my $stop = $tmp[2];
	my $vsinput_l10p = $tmp[3];
	my $vsinput_l2fenr = $tmp[4];
	my $peak_id = $chr.":".$start."-".$stop.":".$str.":".$vsinput_l10p.":".$vsinput_l2fenr;

	push @{$peak_hash{$chr}{$str}{$fi}},$peak_id;
	$peaks2start{$chr}{$str}{$fi}{$peak_id} = $start;
	$peaks2l10p{$chr}{$str}{$fi}{$peak_id} = $vsinput_l10p;


	my $x = int($start / $hashing_value);
	my $y = int( $stop / $hashing_value);

	for my $i ($x..$y) {
	    push @{$read_hash{$chr}{$str}{$fi}{$i}},$peak_id


	}

    }
    close(F);
}
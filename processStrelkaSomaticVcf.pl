use warnings;
use strict;
use Getopt::Long;



my %opt = (
		   "vcf" => undef,
		   "type" => undef,
		   "out" => undef,
		   "DP" => 0);


GetOptions("vcf=s" => \$opt{"vcf"},
		   "type=s" => \$opt{"type"},
		   "out=s" => \$opt{"out"},
		   "dp=i" => \$opt{"DP"});

die "Error: Invalid options\n" if(!defined $opt{"vcf"} || !defined $opt{"type"} || !defined $opt{"out"});


open my $out, ">", $opt{"out"} || die "Error: Cannot open output file.\n";
open my $vcf, "<", $opt{"vcf"} || die "Error: Cannot open input file.\n";

while(<$vcf>){
	chomp;
	if(m/^#/){
		print $out $_, "\n";
		next;
	}
	my @rec = split /\t/;
	next if $rec[6] ne "PASS";
	
	my ($ref, $alt) = ($rec[3], $rec[4]);
	my ($nt, $sgt) = $rec[7] =~ m/NT\=(.+?)\;.+?SGT\=(.+?)\;/;
	
	next if $nt eq 'hom';
	my @sgt = split /->/,$sgt;
	next if $sgt[0] eq $sgt[1];
	

	my($n_gt, $t_gt);	
	if($opt{"type"} =~ m/snv/i){
		my @n_sgt = split //,$sgt[0];
		my @t_sgt = split //,$sgt[1];
		my $ndiff_count = grep{$_ ne $ref}(@n_sgt);
		my $tdiff_count = grep{$_ ne $ref}(@t_sgt);
		
		next if $tdiff_count == 0;
		
		if($ndiff_count == 0){
			$n_gt = '0/0';
		}elsif($ndiff_count == 1){
			$n_gt = '0/1';
		}else{
			$n_gt = '1/1';
		}
		
	
		if($tdiff_count == 0){
			$t_gt = '0/0';
		}elsif($tdiff_count == 1){
			$t_gt = '0/1';
		}else{
			$t_gt = '1/1';
		}		
	}elsif($opt{"type"} =~ m/indel/i){
		if($sgt[0] eq 'het'){
			$n_gt = '0/1';
		}elsif($sgt[0] eq 'ref'){
			$n_gt = '0/0';
		}elsif($sgt[0] eq 'hom'){
			$n_gt = '1/1';
		}
		
		if($sgt[1] eq 'het'){
			$t_gt = '0/1';
		}elsif($sgt[1] eq 'ref'){
			$t_gt = '0/0';
		}elsif($sgt[1] eq 'hom'){
			$t_gt = '1/1';
		}
	}
	
	## Check DP ##
	next if((split /\:/, $rec[9])[0] < $opt{"DP"} || (split /\:/, $rec[10])[0] < $opt{"DP"});
	
	$rec[8] = join(':', ("GT", (split /\:/, $rec[8])));
	$rec[9] = join(':', ($n_gt, (split /\:/, $rec[9])));
	$rec[10] = join(':', ($t_gt, (split /\:/, $rec[10])));
	
	print $out join("\t", @rec), "\n";
}
close $vcf;
close $out;

exit 0;
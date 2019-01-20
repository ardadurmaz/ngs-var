use warnings;
use strict;


my $vcf_file = shift @ARGV;
open my $fh, "<", $vcf_file || die $!;
while(<$fh>){
	chomp;
	if($_=~m/^#/){
		print $_, "\n";
		next;
	}
	my @rec = split /\t/;
	next if $rec[6] ne "PASS" || $rec[5] eq "." || $rec[5] < 10;
	my $flag = 0;
	map{
		my @local_rec = split /:/,$rec[$_];
		$flag = 1 if $local_rec[3] < 8;
	}(9..@rec-1);

	next if $flag == 1;
	print $_, "\n";
}
close $fh;


exit;

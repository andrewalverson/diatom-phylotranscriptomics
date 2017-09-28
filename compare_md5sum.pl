#!/usr/bin/perl

use strict;
use warnings;

my $usage = "\nUsage: $0 provided.md5 post-download.md5\n\n";
$ARGV[1] or die $usage;

my( %original, %download );

# open and store provided md5 file
open( F, $ARGV[0] ) || die "Couldn't open $ARGV[0]: $!\n";
while ( <F>){
  chomp;
  my @a = split /\s+/;
  my @b = split /\//, $a[1];
  $original{pop @b} = $a[0];
}
close F;


# open and store md5 file created after downloading/moving the data
open( F, $ARGV[1] ) || die "Couldn't open $ARGV[1]: $!\n";
while ( <F>){
  chomp;
  my @a = split /\s+/;
  my @b = split /\//, $a[1];
  $download{pop @b} = $a[0];
}
close F;

foreach my $file ( sort keys %original ){
  if( $original{$file} eq $download{$file} ){
    print "$file - OK\n";
  }else{
    print "$file: $original{$file}  ne  $download{$file}\n";
  }
}

exit;

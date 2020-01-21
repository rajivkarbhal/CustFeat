#!/usr/bin/perl
use CustFeat;
print "Enter UniprotKB accession number:\n";
$id = <STDIN>; #Uniprot Accession ID from User
chomp $id;

$seqch=CustFeat::ActSite($id);
print $seqch;

$seqch=CustFeat::Domain($id);
print $seqch;
 
#!/usr/bin/perl

#************** Perl Module for customized extraction of featurer/annotations sequence of UniprotKB entry/ies. **************#
#																															 #
#				Bioinformatics Centre, Savitribai Phule Pune University (Formerly University of Pune)                        #
#                                                                                                                            #        
#****************************************************************************************************************************#


package CustFeat;

use WWW::Mechanize;
$mech = WWW::Mechanize->new();

#$mech->proxy(['http','https'], 'http://username:password@x.x.x.x:xxxx'); Provide Proxy authentication if your network have proxy connection


$emptstrg="";
$dot=".";

#*************************** Subroutine for Sequence Start ***************************
sub seq {  

$id = shift;
chomp $id;

$url = "https://www.uniprot.org/uniprot/?query=id:$id&format=tab&columns=id";
$mech->get($url);
$output = $mech->content();
@data=split ("\n",$output);
$l=@data;

	if ($l==0) 
	{
		print "Please Check Uniprot Accession Number\n";
	}
	else
	{

	# ********** Get Protein Sequence ********** 
	sub Geteq {
	 $sequrl="https://www.uniprot.org/uniprot/$id.fasta";
	 $mech->get($sequrl);
	 $seqout = $mech->content();
	 @sqln=split("\n",$seqout);
		foreach $hhh (@sqln) {
			#print $hhh;
				if ($hhh!~/^\>/) {
					$conc.=$hhh;
			}
		}
		$seq=$conc;
	}

	}
	# ********** Get Protein Sequence Ends********** 


}        

#*************************** Subroutine for Sequence End ***************************



sub ActSite {
open (fw,">cont.txt");
print "Data for Active Site:\n";
$id = shift;
chomp $id;
Geteq($id);

								$feturl="https://www.uniprot.org/uniprot/?query=id:$id&format=tab&columns=feature(ACTIVE%20SITE)";
								$mech->get($feturl);
								$featcontout = $mech->content();
								@flcont=split("\n",$featcontout);
								foreach $hbk (@flcont) {
									chomp $hbk;
									if ($hbk=~/^\S+/) {
										print fw "$hbk\n";
									}
								}

								$out=`findstr /R /N "^" cont.txt | find /C ":"`;

								if ($out>1) {
								 @actdt=split("\n",$featcontout);
								 foreach $acdtln (@actdt) {
									 chomp $acdtln;
									  if ($acdtln=~/^ACT_SITE/g) {										  
											 @splact=split(";",$acdtln);											  
											 foreach $echac (@splact) {												 
												 $echac=~s/^\s+//g;
												 $echac=~s/\s+$//g;
												 if ($echac=~/ACT_SITE\s+(\d+)/) {
													 $stactst=$1;
													 $st=$stactst-1;
													 $temp = substr($seq, $st,1);
													 $out= ">$id | Active Site | $stactst\n$temp\n";
													 print $out,"\n";
												 }
											 }
										 }
								 }
								 }

								else
								{
									print "Annotation Active Site is not present in the entry: $id\n"; 
								}
close (fw);
unlink ("cont.txt");

}

sub AltSeq {
open (fw,">cont.txt");
print "Data for Alternative sequence:\n";
$id = shift;
chomp $id;
Geteq($id);
								$feturl="https://www.uniprot.org/uniprot/?query=id:$id&format=tab&columns=feature(ALTERNATIVE%20SEQUENCE)";
								$mech->get($feturl);
								$featcontout = $mech->content();
								@flcont=split("\n",$featcontout);
								foreach $hbk (@flcont) {
									chomp $hbk;
									if ($hbk=~/^\S+/) {
										print fw "$hbk\n";
									}
								}

								$out=`findstr /R /N "^" cont.txt | find /C ":"`;
								if ($out>1) {
								@actdt=split("\n",$featcontout);
								 foreach $acdtln (@actdt) {
									 chomp $acdtln;									 
									  if ($acdtln=~/^VAR_SEQ/g) {										  
											 @splact=split(";",$acdtln);											  
											 foreach $echac (@splact) {
												 $echac=~s/^\s+//g;
												 $echac=~s/\s+$//g;												 
												 if ($echac=~/VAR_SEQ\s+(\d+)\.\.(\d+)/) {
													 $stactst=$1;
													 $end=$2;
													 $tot=$end-$stactst;
													 $totseq=$tot+1;
													 $st=$stactst-1;
													 $temp = substr($seq, $st,$totseq);
													 $out= ">$id | Alternative sequence | $stactst-$end\n$temp\n";
													 print $out,"\n";
												 }
											 }
										 }
								 }
								 }
								else
								{
									print "Annotation Alternative sequence is not present in the entry: $id\n"; 
								}
close (fw);
unlink ("cont.txt");

}

sub BetaStrand {
open (fw,">cont.txt");
print "Data for Beta Strand:\n";
$id = shift;
chomp $id;
Geteq($id);
								$feturl="https://www.uniprot.org/uniprot/?query=id:$id&format=tab&columns=feature(BETA%20STRAND)";
								$mech->get($feturl);
								$featcontout = $mech->content();
								@flcont=split("\n",$featcontout);
								foreach $hbk (@flcont) {
									chomp $hbk;
									if ($hbk=~/^\S+/) {
										print fw "$hbk\n";
									}
								}

								$out=`findstr /R /N "^" cont.txt | find /C ":"`;
								if ($out>1) {
								@actdt=split("\n",$featcontout);
								 foreach $acdtln (@actdt) {
									 chomp $acdtln;									 
									  if ($acdtln=~/^STRAND/g) {										  
											 @splact=split(";",$acdtln);											  
											 foreach $echac (@splact) {
												 $echac=~s/^\s+//g;
												 $echac=~s/\s+$//g;												 
												 if ($echac=~/STRAND\s+(\d+)\.\.(\d+)/) {
													 $stactst=$1;
													 $end=$2;
													 $tot=$end-$stactst;
													 $totseq=$tot+1;
													 $st=$stactst-1;
													 $temp = substr($seq, $st,$totseq);
													 $out= ">$id | Beta Strand | $stactst-$end\n$temp\n";
													 print $out,"\n";

												 }
											 }
										 }
								 }
								 }
								else
								{
									print "Annotation Beta Strand is not present in the entry: $id\n"; 
								}
close (fw);
unlink ("cont.txt");
}

sub BindSite { 
open (fw,">cont.txt");
print "Data for Binding site:\n";
$id = shift;
chomp $id;
Geteq($id);
								$feturl="https://www.uniprot.org/uniprot/?query=id:$id&format=tab&columns=feature(BINDING%20SITE)";
								$mech->get($feturl);
								$featcontout = $mech->content();								
								@actdt=split("\n",$featcontout);
								@flcont=split("\n",$featcontout);
								foreach $hbk (@flcont) {
									chomp $hbk;
									if ($hbk=~/^\S+/) {
										print fw "$hbk\n";
									}
								}

								$out=`findstr /R /N "^" cont.txt | find /C ":"`;
								if ($out>1) {
								 foreach $acdtln (@actdt) {
									 chomp $acdtln;									 
									  if ($acdtln=~/^BINDING/g) {										  
											 @splact=split(";",$acdtln);											  
											 foreach $echac (@splact) {
												 $echac=~s/^\s+//g;
												 $echac=~s/\s+$//g;												 
												 if ($echac=~/^BINDING\s+(\d+)/) {
													 $stactst=$1;
													 $st=$stactst-1;
													 $temp = substr($seq, $st,1);
													 $out= ">$id | Binding site | $stactst\n$temp\n";
													 print $out,"\n";

												 }
											 }
										 }
								 }
								 }
								else
								{
									print "Annotation Binding site is not present in the entry: $id\n"; 
								}
								close (fw);
unlink ("cont.txt");
}

sub Chain { 
open (fw,">cont.txt");
print "Data for Chain:\n";
$id = shift;
chomp $id;
Geteq($id);
								$feturl="https://www.uniprot.org/uniprot/?query=id:$id&format=tab&columns=feature(CHAIN)";
								$mech->get($feturl);
								$featcontout = $mech->content();
								
								@flcont=split("\n",$featcontout);
								foreach $hbk (@flcont) {
									chomp $hbk;
									if ($hbk=~/^\S+/) {
										print fw "$hbk\n";
									}
								}

								$out=`findstr /R /N "^" cont.txt | find /C ":"`;
								if ($out>1) {
																
								 foreach $acdtln (@actdt) {
									 chomp $acdtln;
									 $acdtln=~s/^\s+//g;
									 $acdtln=~s/\s+$//g;

									  if ($acdtln=~/^CHAIN/g) {
										  chomp $acdtln;
										  
											 @splact=split("\/",$acdtln);
											 foreach $echac (@splact) {												
												 $echac=~s/^\s+//g;
												 $echac=~s/\s+$//g;
												 if ($echac=~/CHAIN\s+(\d+)\.\.(\d+)/) {
													 $stactst=$1;
													 $st=$stactst-1;
													 $end=$2;
													 $tot=$end-$stactst;
													 $totseq=$tot+1;
													 $st=$stactst-1;
													 $temp = substr($seq, $st,$totseq);
													 $out= ">$id | Chain | $stactst-$end\n$temp\n";
													 print $out,"\n";
												 }
											 }
										 }									
								 
								}
								}
								else
								{
									print "Annotation Chain is not present in the entry: $id\n"; 
								}
close (fw);
unlink ("cont.txt");								
}


sub CoilCoil { 
open (fw,">cont.txt");
print "Data for Coiled coil:\n";
$id = shift;
chomp $id;
Geteq($id);
								$feturl="https://www.uniprot.org/uniprot/?query=id:$id&format=tab&columns=feature(COILED%20COIL)";
								$mech->get($feturl);
								$featcontout = $mech->content();
								@flcont=split("\n",$featcontout);
								foreach $hbk (@flcont) {
									chomp $hbk;
									if ($hbk=~/^\S+/) {
										print fw "$hbk\n";
									}
								}

								$out=`findstr /R /N "^" cont.txt | find /C ":"`;
								if ($out>1) {
								@actdt=split("\n",$featcontout);
								 foreach $acdtln (@actdt) {
									 chomp $acdtln;									 
									  if ($acdtln=~/^COILED/g) {										  
											 @splact=split(";",$acdtln);											  
											 foreach $echac (@splact) {
												 $echac=~s/^\s+//g;
												 $echac=~s/\s+$//g;												 
												 if ($echac=~/COILED\s+(\d+)\.\.(\d+)/) {
													 $stactst=$1;
													 $st=$stactst-1;
													 $end=$2;
													 $tot=$end-$stactst;
													 $totseq=$tot+1;
													 $st=$stactst-1;
													 $temp = substr($seq, $st,$totseq);
													 $out= ">$id | Coiled coil | $stactst-$end\n$temp\n";
													 print $out,"\n";
												 }
											 }
										 }
								 }
								 }
								else
								{
									print "Annotation Domain is not present in the entry: $id\n"; 
								}
close (fw);
unlink ("cont.txt");
}


sub CompBias {
open (fw,">cont.txt");
print "Data for Compositional bias:\n";
$id = shift;
chomp $id;
Geteq($id);
								$feturl="https://www.uniprot.org/uniprot/?query=id:$id&format=tab&columns=feature(COMPOSITIONAL%20BIAS)";
								$mech->get($feturl);
								$featcontout = $mech->content();
								@flcont=split("\n",$featcontout);
								foreach $hbk (@flcont) {
									chomp $hbk;
									if ($hbk=~/^\S+/) {
										print fw "$hbk\n";
									}
								}

								$out=`findstr /R /N "^" cont.txt | find /C ":"`;
								if ($out>1) {
								@actdt=split("\n",$featcontout);
								 foreach $acdtln (@actdt) {
									 chomp $acdtln;									 
									  if ($acdtln=~/^COMPBIAS/g) {										  
											 @splact=split(";",$acdtln);											  
											 foreach $echac (@splact) {
												 $echac=~s/^\s+//g;
												 $echac=~s/\s+$//g;												 
												 if ($echac=~/COMPBIAS\s+(\d+)\.\.(\d+)/) {
													 $stactst=$1;
													 $st=$stactst-1;
													 $end=$2;
													 $tot=$end-$stactst;
													 $totseq=$tot+1;
													 $st=$stactst-1;
													 $temp = substr($seq, $st,$totseq);
													 $out= ">$id | Compositional bias | $stactst-$end\n$temp\n";
													 print $out,"\n";
												 }
											 }
										 }
								 }
								 }
								else
								{
									print "Annotation Compositional bias is not present in the entry: $id\n"; 
								}
close (fw);
unlink ("cont.txt");
}


sub Crosslink { 
open (fw,">cont.txt");
print "Data for Cross-link:\n";
$id = shift;
chomp $id;
Geteq($id);
								$feturl="https://www.uniprot.org/uniprot/?query=id:$id&format=tab&columns=feature(CROSS%20LINK)";
								$mech->get($feturl);
								$featcontout = $mech->content();
								@flcont=split("\n",$featcontout);
								foreach $hbk (@flcont) {
									chomp $hbk;
									if ($hbk=~/^\S+/) {
										print fw "$hbk\n";
									}
								}

								$out=`findstr /R /N "^" cont.txt | find /C ":"`;
								if ($out>1) {
								@actdt=split("\n",$featcontout);
								 foreach $acdtln (@actdt) {
									 chomp $acdtln;									 
									  if ($acdtln=~/^CROSSLNK/g) {										  
											 @splact=split(";",$acdtln);											  
											 foreach $echac (@splact) {
												 $echac=~s/^\s+//g;
												 $echac=~s/\s+$//g;												 
												 if ($echac=~/CROSSLNK\s+(\d+)/) {
													 $stactst=$1;
													 $st=$stactst-1;
													 $temp = substr($seq, $st,1);
													 $out= ">$id | Cross-link | $stactst\n$temp\n";
													 print $out,"\n";
												 }
											 }
										 }
								 }
								 }
								else
								{
									print "Annotation Cross-link is not present in the entry: $id\n"; 
								}
close (fw);
unlink ("cont.txt");
}


sub DisulBond { 
open (fw,">cont.txt");
print "Data for Disulfide bond:\n";
$id = shift;
chomp $id;
Geteq($id);
								$feturl="https://www.uniprot.org/uniprot/?query=id:$id&format=tab&columns=feature(DISULFIDE BOND)";
								$mech->get($feturl);
								$featcontout = $mech->content();
								@flcont=split("\n",$featcontout);
								foreach $hbk (@flcont) {
									chomp $hbk;
									if ($hbk=~/^\S+/) {
										print fw "$hbk\n";
									}
								}

								$out=`findstr /R /N "^" cont.txt | find /C ":"`;
								if ($out>1) {
								@actdt=split("\n",$featcontout);
								 foreach $acdtln (@actdt) {
									 chomp $acdtln;									 
									  if ($acdtln=~/^DISULFID/g) {										  
											 @splact=split(";",$acdtln);											  
											 foreach $echac (@splact) {
												 $echac=~s/^\s+//g;
												 $echac=~s/\s+$//g;
												 if ($echac=~/DISULFID\s+(\d+)\.\.(\d+)/) {
													 $stactst=$1;
													 $st=$stactst-1;
													 $end=$2;
													 $ed=$end-1;
													 @dsbd=split("",$seq);
													 $out= ">$id | Disulfide bond | $stactst-$end\n$dsbd[$st] <-> $dsbd[$ed]\n";
													 print $out,"\n";
												 }
											 }
										 }
								 }
								 }
								else
								{
									print "Annotation Disulfide bond is not present in the entry: $id\n"; 
								}
close (fw);
unlink ("cont.txt");
								 
}


sub Domain {
open (fw,">cont.txt");
print "Data for Domain:\n";
$id = shift;
chomp $id;
Geteq($id);
								$feturl="https://www.uniprot.org/uniprot/?query=id:$id&format=tab&columns=feature(DOMAIN EXTENT)";
								$mech->get($feturl);
								$featcontout = $mech->content();
								@flcont=split("\n",$featcontout);
								foreach $hbk (@flcont) {
									chomp $hbk;
									if ($hbk=~/^\S+/) {
										print fw "$hbk\n";
									}
								}

								$out=`findstr /R /N "^" cont.txt | find /C ":"`;
								if ($out>1) {
								@actdt=split("\n",$featcontout);
								 foreach $acdtln (@actdt) {
									 chomp $acdtln;	
									 $acdtln=~s/\s+DOMAIN/__DOMAIN/g;
									  if ($acdtln=~/^DOMAIN/g) {
										  
											 @splact=split("__",$acdtln);
											 										 
											 foreach $echac (@splact) {	
												 #print $echac,"\n";
												 $echac=~s/^\s+//g;
												 $echac=~s/\s+$//g;	
												 
												 if ($echac=~/^DOMAIN\s+(\d+)\.\.(\d+);\s+(.+)/) {
													 $stactst=$1;
													 $st=$stactst-1;
													 $end=$2;
													 $dnm=$3;
													 $dnm=~s/\/note\=\"//g;
													 $dnm=~s/\"\;//g;
													 $dnm=~s/\s+\/evidence=.+//g;
													 $tot=$end-$stactst;
													 $totseq=$tot+1;
													 $st=$stactst-1;
													 $temp = substr($seq, $st,$totseq);
													 $out= ">$id | Domain: $dnm | $stactst-$end\n$temp\n";
													 print $out,"\n";
												 }
											 }
																			  
										 }
								 }
								 }
								else
								{
									print "Annotation Domain is not present in the entry: $id\n"; 
								}
close (fw);
unlink ("cont.txt");
}

sub Glyco {
open (fw,">cont.txt");
print "Data for Glycosylation:\n";
$id = shift;
chomp $id;
Geteq($id);
								$feturl="https://www.uniprot.org/uniprot/?query=id:$id&format=tab&columns=feature(GLYCOSYLATION)";
								$mech->get($feturl);
								$featcontout = $mech->content();
								@flcont=split("\n",$featcontout);
								foreach $hbk (@flcont) {
									chomp $hbk;
									if ($hbk=~/^\S+/) {
										print fw "$hbk\n";
									}
								}

								$out=`findstr /R /N "^" cont.txt | find /C ":"`;
								if ($out>1) {
								@actdt=split("\n",$featcontout);
								 foreach $acdtln (@actdt) {
									 chomp $acdtln;	
									 $acdtln=~s/\s+CARBOHYD/__CARBOHYD/g;
									  if ($acdtln=~/^CARBOHYD/g) {										  
											 @splact=split("__",$acdtln);											  
											 foreach $echac (@splact) {
												 $echac=~s/^\s+//g;
												 $echac=~s/\s+$//g;
												 if ($echac=~/CARBOHYD\s+(\d+)\;(.+)/) {
													 $stactst=$1;
													 $st=$stactst-1;
													 $dnm=$2;
													 $dnm=~s/\/note\=\"//g;
													 $dnm=~s/\"\;//g;
													 $dnm=~s/\s+\/evidence=.+//g;
													 $ed=$end-1;
													 @dsbd=split("",$seq);
													 $temp = substr($seq, $st,1);
													 $out= ">$id | Glycosylation: $dnm | $stactst\n$temp\n";
													 print $out,"\n";
												 }
											 }
										 }
								 }
								 }
								else
								{
									print "Annotation Glycosylation is not present in the entry: $id\n"; 
								}
close (fw);
unlink ("cont.txt");
}

sub Initmeth { 
open (fw,">cont.txt");
print "Data for Initiator methionine:\n";
$id = shift;
chomp $id;
Geteq($id);
								$feturl="https://www.uniprot.org/uniprot/?query=id:$id&format=tab&columns=feature(INITIATOR%20METHIONINE)";
								$mech->get($feturl);
								$featcontout = $mech->content();
								@flcont=split("\n",$featcontout);
								foreach $hbk (@flcont) {
									chomp $hbk;
									if ($hbk=~/^\S+/) {
										print fw "$hbk\n";
									}
								}

								$out=`findstr /R /N "^" cont.txt | find /C ":"`;
								if ($out>1) {
								@actdt=split("\n",$featcontout);
								 foreach $acdtln (@actdt) {
									 chomp $acdtln;									 
									  if ($acdtln=~/^INIT_MET/g) {										  
											 @splact=split("\;",$acdtln);											  
											 foreach $echac (@splact) {
												 $echac=~s/^\s+//g;
												 $echac=~s/\s+$//g;
												# print $echac,"\n";
												 if ($echac=~/INIT_MET\s+(\d+)\s+(\d+)\s+(.+)./) {
													 $stactst=$1;
													 $st=$stactst-1;
													 $dnm=$3;
													 $dnm=~s/ \{.+\}//g;
													 $temp = substr($seq, $st,1);
													 $out= ">$id | Initiator methionine: $dnm | $stactst\n$temp\n";
													 print $out,"\n";
												 }
											 }
										 }
								 }
								 }
								else
								{
									print "Annotation Initiator methionine is not present in the entry: $id\n"; 
								}
close (fw);
unlink ("cont.txt");
}

sub Intramemb {
open (fw,">cont.txt");	
print "Data for Intramembrane:\n";
$id = shift;
chomp $id;
Geteq($id);
								$feturl="https://www.uniprot.org/uniprot/?query=id:$id&format=tab&columns=feature(INTRAMEMBRANE)";
								$mech->get($feturl);
								$featcontout = $mech->content();
								@flcont=split("\n",$featcontout);
								foreach $hbk (@flcont) {
									chomp $hbk;
									if ($hbk=~/^\S+/) {
										print fw "$hbk\n";
									}
								}

								$out=`findstr /R /N "^" cont.txt | find /C ":"`;
								if ($out>1) {
								@actdt=split("\n",$featcontout);
								 foreach $acdtln (@actdt) {
									 chomp $acdtln;									 
									  if ($acdtln=~/^INTRAMEM/g) {										  
											 @splact=split("\;",$acdtln);											  
											 foreach $echac (@splact) {
												 $echac=~s/^\s+//g;
												 $echac=~s/\s+$//g;	
												 if ($echac=~/INTRAMEM\s+(\d+)\.\.(\d+)/) {
													 $stactst=$1;
													 $st=$stactst-1;
													 $end=$2;
													 $dnm=$3;
													 $dnm=~s/ \{.+\}//g;
													 $tot=$end-$stactst;
													 $totseq=$tot+1;
													 $temp = substr($seq, $st,$totseq);
													 $out= ">$id | Intramembrane | $stactst-$end\n$temp\n";
													 print $out,"\n";
												 }
											 }
										 }
								 }
								 }
								else
								{
									print "Annotation Domain is not present in the entry: $id\n"; 
								}
close (fw);
unlink ("cont.txt");
}

sub Lipidation {
open (fw,">cont.txt");
print "Data for Lipidation:\n";
$id = shift;
chomp $id;
Geteq($id);
								$feturl="https://www.uniprot.org/uniprot/?query=id:$id&format=tab&columns=feature(LIPIDATION)";
								$mech->get($feturl);
								$featcontout = $mech->content();
								@flcont=split("\n",$featcontout);
								foreach $hbk (@flcont) {
									chomp $hbk;
									if ($hbk=~/^\S+/) {
										print fw "$hbk\n";
									}
								}

								$out=`findstr /R /N "^" cont.txt | find /C ":"`;
								if ($out>1) {
								@actdt=split("\n",$featcontout);
								 foreach $acdtln (@actdt) {
									 chomp $acdtln;									 
									  if ($acdtln=~/^LIPID/g) {
										   @splact=split("\;",$acdtln);											  
											 foreach $echac (@splact) {
												 $echac=~s/^\s+//g;
												 $echac=~s/\s+$//g;
												 if ($echac=~/LIPID\s+(\d+)/) {
													 $stactst=$1;
													 $st=$stactst-1;
													 $temp = substr($seq, $st,1);
													 $out= ">$id | Lipidation | $stactst\n$temp\n";
													 print $out,"\n";
												 }
											 }
										 }
								 }
								 }
								else
								{
									print "Annotation Domain is not present in the entry: $id\n"; 
								}
close (fw);
unlink ("cont.txt");
}

sub MetBind { 
open (fw,">cont.txt");
print "Data for Metal binding:\n";
$id = shift;
chomp $id;
Geteq($id);
								$feturl="https://www.uniprot.org/uniprot/?query=id:$id&format=tab&columns=feature(METAL%20BINDING)";
								$mech->get($feturl);
								$featcontout = $mech->content();
								@flcont=split("\n",$featcontout);
								foreach $hbk (@flcont) {
									chomp $hbk;
									if ($hbk=~/^\S+/) {
										print fw "$hbk\n";
									}
								}

								$out=`findstr /R /N "^" cont.txt | find /C ":"`;
								if ($out>1) {
								@actdt=split("\n",$featcontout);
								 foreach $acdtln (@actdt) {
									 chomp $acdtln;									 
									  if ($acdtln=~/^METAL/g) {										  
											 @splact=split("\;",$acdtln);											  
											 foreach $echac (@splact) {
												 $echac=~s/^\s+//g;
												 $echac=~s/\s+$//g;
												 if ($echac=~/METAL\s+(\d+)/) {
													 $stactst=$1;
													 $st=$stactst-1;
													 $dnm=$3;
													 $dnm=~s/ \{.+\}//g;
													 $temp = substr($seq, $st,1);
													 $out= ">$id | Metal binding | $stactst\n$temp\n";
													 print $out,"\n";
												 }
											 }
										 }
								 }
								 }
								else
								{
									print "Annotation Metal binding is not present in the entry: $id\n"; 
								}
close (fw);
unlink ("cont.txt");
}

sub DNABind { 
open (fw,">cont.txt");
print "Data for DNA binding:\n";
$id = shift;
chomp $id;
Geteq($id);
								$feturl="https://www.uniprot.org/uniprot/?query=id:$id&format=tab&columns=feature(DNA BINDING)";
								$mech->get($feturl);
								$featcontout = $mech->content();
								@flcont=split("\n",$featcontout);
								foreach $hbk (@flcont) {
									chomp $hbk;
									if ($hbk=~/^\S+/) {
										print fw "$hbk\n";
									}
								}

								$out=`findstr /R /N "^" cont.txt | find /C ":"`;
								if ($out>1) {
								@actdt=split("\n",$featcontout);
								 foreach $acdtln (@actdt) {
									 chomp $acdtln;									 
									  if ($acdtln=~/^DNA_BIND/g) {
										  @splact=split("\;",$acdtln);											  
											 foreach $echac (@splact) {
												 $echac=~s/^\s+//g;
												 $echac=~s/\s+$//g;												 
												 if ($echac=~/DNA_BIND\s+(\d+)\.\.(\d+)/) {
													 $stactst=$1;
													 $st=$stactst-1;
													 $end=$2;
													 $tot=$end-$stactst;
													 $totseq=$tot+1;
													 $st=$stactst-1;
													 $temp = substr($seq, $st,$totseq);
													 $out= ">$id | DNA binding| $stactst-$end\n$temp\n";
													 print $out,"\n";
												 }
											 }
										 }
								 }
								 }
								else
								{
									print "Annotation DNA binding is not present in the entry: $id\n"; 
								}
close (fw);
unlink ("cont.txt");
}

sub Helix {
open (fw,">cont.txt");	
print "Data for Helix\n";
$id = shift;
chomp $id;
Geteq($id);
								$feturl="https://www.uniprot.org/uniprot/?query=id:$id&format=tab&columns=feature(HELIX)";
								$mech->get($feturl);
								$featcontout = $mech->content();
								@flcont=split("\n",$featcontout);
								foreach $hbk (@flcont) {
									chomp $hbk;
									if ($hbk=~/^\S+/) {
										print fw "$hbk\n";
									}
								}

								$out=`findstr /R /N "^" cont.txt | find /C ":"`;
								if ($out>1) {
								@actdt=split("\n",$featcontout);
								 foreach $acdtln (@actdt) {
									 chomp $acdtln;									 
									  if ($acdtln=~/^HELIX/g) {										  
											 @splact=split("\;",$acdtln);											  
											 foreach $echac (@splact) {
												 $echac=~s/^\s+//g;
												 $echac=~s/\s+$//g;												 
												 if ($echac=~/HELIX\s+(\d+)\.\.(\d+)/) {
													 $stactst=$1;
													 $st=$stactst-1;
													 $end=$2;
													 $tot=$end-$stactst;
													 $totseq=$tot+1;
													 $temp = substr($seq, $st,$totseq);
													 $out= ">$id | Helix | $stactst-$end\n$temp\n";
													 print $out,"\n";
												 }
											 }
										 }
								 }
								 }
								else
								{
									print "Annotation Helix is not present in the entry: $id\n"; 
								}
close (fw);
unlink ("cont.txt");
}

sub ModRes {
open (fw,">cont.txt");	
print "Data for Modified residue:\n";
$id = shift;
chomp $id;
Geteq($id);
								$feturl="https://www.uniprot.org/uniprot/?query=id:$id&format=tab&columns=feature(MODIFIED%20RESIDUE)";
								$mech->get($feturl);
								$featcontout = $mech->content();
								@flcont=split("\n",$featcontout);
								foreach $hbk (@flcont) {
									chomp $hbk;
									if ($hbk=~/^\S+/) {
										print fw "$hbk\n";
									}
								}

								$out=`findstr /R /N "^" cont.txt | find /C ":"`;
								if ($out>1) {
								@actdt=split("\n",$featcontout);
								 foreach $acdtln (@actdt) {
									 chomp $acdtln;	
									  $acdtln=~s/\s+MOD_RES/__MOD_RES/g;
									  if ($acdtln=~/^MOD_RES/g) {										  
											 @splact=split("__",$acdtln);											  
											 foreach $echac (@splact) {
												 $echac=~s/^\s+//g;
												 $echac=~s/\s+$//g;
												 if ($echac=~/MOD_RES\s+(\d+);\s+(.+)/) {
													 $stactst=$1;
													 $st=$stactst-1;
													 $dnm=$2;
													 $dnm=~s/\/note\=\"//g;
													 $dnm=~s/\"\;//g;
													 $dnm=~s/\s+\/evidence=.+//g;
													 $ed=$end-1;
													 @dsbd=split("",$seq);

													 $out= ">$id | Modified residue: $dnm | $stactst\n$dsbd[$st]\n";
													 print $out,"\n";
												 }
											 }
										 }
								 }
								 }
								else
								{
									print "Annotation Modified residue is not present in the entry: $id\n"; 
								}
close (fw);
unlink ("cont.txt");
}


sub Motif { 
open (fw,">cont.txt");
print "Data for Motif:\n";
$id = shift;
chomp $id;
Geteq($id);
								$feturl="https://www.uniprot.org/uniprot/?query=id:$id&format=tab&columns=feature(MOTIF)";
								$mech->get($feturl);
								$featcontout = $mech->content();
								@flcont=split("\n",$featcontout);
								foreach $hbk (@flcont) {
									chomp $hbk;
									if ($hbk=~/^\S+/) {
										print fw "$hbk\n";
									}
								}

								$out=`findstr /R /N "^" cont.txt | find /C ":"`;
								if ($out>1) {
								@actdt=split("\n",$featcontout);
								 foreach $acdtln (@actdt) {
									 chomp $acdtln;
									 $acdtln=~s/\s+MOTIF/__MOTIF/g;
									  if ($acdtln=~/^MOTIF/g) {										  
											 @splact=split("__",$acdtln);											  
											 foreach $echac (@splact) {
												 $echac=~s/^\s+//g;
												 $echac=~s/\s+$//g;												 
												 if ($echac=~/MOTIF\s+(\d+)\.\.(\d+)\;\s+(.+)/) {
													 $stactst=$1;
													 $st=$stactst-1;
													 $end=$2;
													 $dnm=$3;
													 $dnm=~s/\/note\=\"//g;
													 $dnm=~s/\"\;//g;
													 $dnm=~s/\s+\/evidence=.+//g;
													 $tot=$end-$stactst;
													 $totseq=$tot+1;
													 $temp = substr($seq, $st,$totseq);
													 $out= ">$id | Motif: $dnm | $stactst-$end\n$temp\n";
													 print $out,"\n";
												 }
											 }
										 }
								 }
								 }
								else
								{
									print "Annotation Motif is not present in the entry: $id\n"; 
								}
								close (fw);
unlink ("cont.txt");
}


sub NatVar { 
open (fw,">cont.txt");
print "Data for Natural variant:\n";
$id = shift;
chomp $id;
Geteq($id);
								$feturl="https://www.uniprot.org/uniprot/?query=id:$id&format=tab&columns=feature(NATURAL VARIANT)";
								$mech->get($feturl);
								 $featcontout = $mech->content();
								 @flcont=split("\n",$featcontout);
								foreach $hbk (@flcont) {
									chomp $hbk;
									if ($hbk=~/^\S+/) {
										print fw "$hbk\n";
									}
								}

								$out=`findstr /R /N "^" cont.txt | find /C ":"`;
								if ($out>1) {
								 @actdt=split("\n",$featcontout);
								 foreach $acdtln (@actdt) {
									 chomp $acdtln;									 
									  if ($acdtln=~/^VARIANT/g) {
											 @splact=split("\;",$acdtln);											  
											 foreach $echac (@splact) {												 
												 $echac=~s/^\s+//g;
												 $echac=~s/\s+$//g;												 
												 if ($echac=~/VARIANT\s+(\d+)/) {													 
													 $stactst=$1;													 
													 $st=$stactst-1;
													 $temp = substr($seq, $st,1);
													 $out=">$id | Natural variant | $stactst\n$temp\n";
													 print $out,"\n";
												 }
											 }
										 }
								 }
								 }
								else
								{
									print "Annotation Natural variant is not present in the entry: $id\n"; 
								}
close (fw);
unlink ("cont.txt");
}


sub Nonadres { 
open (fw,">cont.txt");
print "Data for Non-adjacent residues:\n";
$id = shift;
chomp $id;
Geteq($id);
								$feturl="https://www.uniprot.org/uniprot/?query=id:$id&format=tab&columns=feature(NON%20ADJACENT%20RESIDUES)";
								$mech->get($feturl);
								 $featcontout = $mech->content();
								 @flcont=split("\n",$featcontout);
								foreach $hbk (@flcont) {
									chomp $hbk;
									if ($hbk=~/^\S+/) {
										print fw "$hbk\n";
									}
								}

								$out=`findstr /R /N "^" cont.txt | find /C ":"`;
								if ($out>1) {
								 @actdt=split("\n",$featcontout);
								 foreach $acdtln (@actdt) {
									 chomp $acdtln;									 
									  if ($acdtln=~/^NON_CONS/g) {										  
											 @splact=split("\;",$acdtln);											  
											 foreach $echac (@splact) {
												 $echac=~s/^\s+//g;
												 $echac=~s/\s+$//g;												 
												 if ($echac=~/NON_CONS\s+(\d+)\.\.(\d+)/) {
													 $stactst=$1;
													 $st=$stactst-1;
													 $end=$2;
													 $dnm=$3;
													 $dnm=~s/ \{.+\}//g;
													 $tot=$end-$stactst;
													 $totseq=$tot+1;
													 $temp = substr($seq, $st,$totseq);
													 $out= ">$id | Non-adjacent residues | $stactst-$end\n$temp\n";
													 print $out,"\n";
												 }
											 }
										 }
								 }
								 }
								else
								{
									print "Annotation Non-adjacent residues is not present in the entry: $id\n"; 
								}
								close (fw);
unlink ("cont.txt");
}


sub Nonstdres { 
open (fw,">cont.txt");
print "Data for Non-standard residue:\n";
$id = shift;
chomp $id;
Geteq($id);
								$feturl="https://www.uniprot.org/uniprot/?query=id:$id&format=tab&columns=feature(NON STANDARD RESIDUE)";
								$mech->get($feturl);
								 $featcontout = $mech->content();
								 @flcont=split("\n",$featcontout);
								foreach $hbk (@flcont) {
									chomp $hbk;
									if ($hbk=~/^\S+/) {
										print fw "$hbk\n";
									}
								}

								$out=`findstr /R /N "^" cont.txt | find /C ":"`;
								if ($out>1) {
								 @actdt=split("\n",$featcontout);
								 foreach $acdtln (@actdt) {
									 chomp $acdtln;									 
									  if ($acdtln=~/^NON_STD/g) {
										  $acdtln=~s/\s+NON_STD/__NON_STD/g;
											 @splact=split("__",$acdtln);											  
											 foreach $echac (@splact) {
												 $echac=~s/^\s+//g;
												 $echac=~s/\s+$//g;
												 if ($echac=~/NON_STD\s+(\d+)\;\s+\/note\=\"(.+)\"\;\s+\/.+/) {
													 $stactst=$1;
													 $st=$stactst-1;
													 $end=$1;
													 $dnm=$2;
													 #$dnm=~s/ \{.+\}//g;
													 $ed=$end-1;
													 @dsbd=split("",$seq);
													 $out= ">$id | Non-standard residue: $dnm | $stactst\n$dsbd[$st]\n";
													 print $out,"\n";
												 }
											 }
										 }
								 }
								 }
								else
								{
									print "Annotation Non-standard residue is not present in the entry: $id\n"; 
								}
close (fw);
unlink ("cont.txt");
}


sub Nonterres { 
print "Data for Non-terminal residue:\n";
$id = shift;
chomp $id;
Geteq($id);
                                $feturl="https://www.uniprot.org/uniprot/?query=id:$id&format=tab&columns=feature(NON%20TERMINAL%20RESIDUE)";
								$mech->get($feturl);
								 $featcontout = $mech->content();								
								 @actdt=split("\n",$featcontout);
								 foreach $acdtln (@actdt) {
									 chomp $acdtln;									 
									  if ($acdtln=~/^NON_TER/g) {										  
											 @splact=split("\;",$acdtln);											  
											 foreach $echac (@splact) {
												 $echac=~s/^\s+//g;
												 $echac=~s/\s+$//g;
												 if ($echac=~/NON_TER\s+(\d+)/) {
													 $stactst=$1;
													 $st=$stactst-1;
													 $end=$1;
													 $ed=$end-1;
													 @dsbd=split("",$seq);
													 $out= ">$id | Non-terminal residue | $stactst\n$dsbd[$st]\n";
													 print $out,"\n";
												 }
											 }
										 }
								 }
}


sub Nubind { 
print "Data for Nucleotide binding:\n";
$id = shift;
chomp $id;
Geteq($id);
								$feturl="https://www.uniprot.org/uniprot/?query=id:$id&format=tab&columns=feature(NP%20BIND)";
								$mech->get($feturl);
								$featcontout = $mech->content();								
								 @actdt=split("\n",$featcontout);
								 foreach $acdtln (@actdt) {
									 chomp $acdtln;									 
									  if ($acdtln=~/^NP_BIND/g) {
										  $acdtln=~s/\s+NP_BIND/__NP_BIND/g;
											 @splact=split("__",$acdtln);											  
											 foreach $echac (@splact) {
												 $echac=~s/^\s+//g;
												 $echac=~s/\s+$//g;												 
												 if ($echac=~/NP_BIND\s+(\d+)\.\.(\d+)\;\s+\/note\=\"(.+)\";.+/) {
													 $stactst=$1;
													 $st=$stactst-1;
													 $end=$2;
													 $dnm=$3;
													 $tot=$end-$stactst;
													 $totseq=$tot+1;
													 $temp = substr($seq, $st,$totseq);
													 $out= ">$id | Nucleotide binding: $dnm | $stactst-$end\n$temp\n";
													 print $out,"\n";
												 }
											 }
										 }
								 }
}


sub Peptide {
open (fw,">cont.txt");
print "Data for Peptide:\n";
$id = shift;
chomp $id;
Geteq($id);
								$feturl="https://www.uniprot.org/uniprot/?query=id:$id&format=tab&columns=feature(PEPTIDE)";
								$mech->get($feturl);
								 $featcontout = $mech->content();
								 @flcont=split("\n",$featcontout);
								foreach $hbk (@flcont) {
									chomp $hbk;
									if ($hbk=~/^\S+/) {
										print fw "$hbk\n";
									}
								}

								$out=`findstr /R /N "^" cont.txt | find /C ":"`;
								if ($out>1) {
								 @actdt=split("\n",$featcontout);
								 foreach $acdtln (@actdt) {
									 chomp $acdtln;									 
									  if ($acdtln=~/^PEPTIDE/g) {
										  $acdtln=~s/\s+PEPTIDE/__PEPTIDE/g;
											 @splact=split("__",$acdtln);											  
											 foreach $echac (@splact) {
												 $echac=~s/^\s+//g;
												 $echac=~s/\s+$//g;												 
												 if ($echac=~/PEPTIDE\s+(\d+)\.\.(\d+)\;\s+\/note\=\"(.+)\"\;\s+/) {
													 $stactst=$1;
													 $st=$stactst-1;
													 $end=$2;
													 $dnm=$3;
													 $dnm=~s/\"\;\s+\/evidence.+//;
													 $tot=$end-$stactst;
													 $totseq=$tot+1;
													 $temp = substr($seq, $st,$totseq);
													 $out= ">$id | Peptide: $dnm | $stactst-$end\n$temp\n";
													 print $out,"\n";
												 }
											 }
										 }
								 }
								 }
								else
								{
									print "Annotation Peptide is not present in the entry: $id\n"; 
								}
								close (fw);
unlink ("cont.txt");
}


sub Propeptide { 
open (fw,">cont.txt");
print "Data for Propeptide:\n";
$id = shift;
chomp $id;
Geteq($id);
								$feturl="https://www.uniprot.org/uniprot/?query=id:$id&format=tab&columns=feature(PROPEPTIDE)";
								$mech->get($feturl);
								 $featcontout = $mech->content();
								 @flcont=split("\n",$featcontout);
								foreach $hbk (@flcont) {
									chomp $hbk;
									if ($hbk=~/^\S+/) {
										print fw "$hbk\n";
									}
								}

								$out=`findstr /R /N "^" cont.txt | find /C ":"`;
								if ($out>1) {
								 @actdt=split("\n",$featcontout);
								 foreach $acdtln (@actdt) {
									 chomp $acdtln;									 
									  if ($acdtln=~/^PROPEP/g) {
										  $acdtln=~s/\s+PROPEP/__PROPEP/g;
											 @splact=split("__",$acdtln);											  
											 foreach $echac (@splact) {
												 $echac=~s/^\s+//g;
												 $echac=~s/\s+$//g;												 
												 if ($echac=~/PROPEP\s+(\d+)\.\.(\d+)/) {
													 $stactst=$1;
													 $st=$stactst-1;
													 $end=$2;
													 $tot=$end-$stactst;
													 $totseq=$tot+1;
													 $temp = substr($seq, $st,$totseq);
													 $out= ">$id | Propeptide | $stactst-$end\n$temp\n";
													 print $out,"\n";
												 }
											 }
										 }
								 }
								 }
								else
								{
									print "Annotation Propeptide is not present in the entry: $id\n"; 
								}
								close (fw);
unlink ("cont.txt");
}


sub Region {
open (fw,">cont.txt");
print "Data for Region:\n";
$id = shift;
chomp $id;
Geteq($id);
								$feturl="https://www.uniprot.org/uniprot/?query=id:$id&format=tab&columns=feature(REGION)";
								$mech->get($feturl);
								 $featcontout = $mech->content();
								 @flcont=split("\n",$featcontout);
								foreach $hbk (@flcont) {
									chomp $hbk;
									if ($hbk=~/^\S+/) {
										print fw "$hbk\n";
									}
								}

								$out=`findstr /R /N "^" cont.txt | find /C ":"`;
								if ($out>1) {
								 @actdt=split("\n",$featcontout);
								 foreach $acdtln (@actdt) {
									 chomp $acdtln;									 
									  if ($acdtln=~/^REGION/g) {
										  $acdtln=~s/\s+REGION/__REGION/g;
											 @splact=split("__",$acdtln);											  
											 foreach $echac (@splact) {
												 $echac=~s/^\s+//g;
												 $echac=~s/\s+$//g;												 
												 if ($echac=~/REGION\s+(\d+)\.\.(\d+)\;\s+\/note\=\"(.+)\"\;/) {
													 $stactst=$1;
													 $st=$stactst-1;
													 $end=$2;
													 $dnm=$3;
													 $dnm=~s/\";\s+\/.+//;
													 $tot=$end-$stactst;
													 $totseq=$tot+1;
													 $temp = substr($seq, $st,$totseq);
													 $out= ">$id | Region: $dnm | $stactst-$end\n$temp\n";
													 print $out,"\n";
												 }
											 }
										 }
								 }
								 }
								else
								{
									print "Annotation Region is not present in the entry: $id\n"; 
								}
								close (fw);
unlink ("cont.txt");
}


sub Repeat { 
open (fw,">cont.txt");
print "Data for Repeat:\n";
$id = shift;
chomp $id;
Geteq($id);
								$feturl="https://www.uniprot.org/uniprot/?query=id:$id&format=tab&columns=feature(REPEAT)";
								$mech->get($feturl);
								 $featcontout = $mech->content();
								 @flcont=split("\n",$featcontout);
								foreach $hbk (@flcont) {
									chomp $hbk;
									if ($hbk=~/^\S+/) {
										print fw "$hbk\n";
									}
								}

								$out=`findstr /R /N "^" cont.txt | find /C ":"`;
								if ($out>1) {
								 @actdt=split("\n",$featcontout);
								 foreach $acdtln (@actdt) {
									 chomp $acdtln;									 
									  if ($acdtln=~/^REPEAT/g) {
										  $acdtln=~s/\s+REPEAT/__REPEAT/g;
											 @splact=split("__",$acdtln);											  
											 foreach $echac (@splact) {
												 $echac=~s/^\s+//g;
												 $echac=~s/\s+$//g;												 
												 if ($echac=~/REPEAT\s+(\d+)\.\.(\d+)\;\s+\/note\=\"(.+)\"\;/) {
													 $stactst=$1;
													 $st=$stactst-1;
													 $end=$2;
													 $dnm=$3;
													 $tot=$end-$stactst;
													 $totseq=$tot+1;
													 $temp = substr($seq, $st,$totseq);
													 $out= ">$id | Repeat: $dnm | $stactst-$end\n$temp\n";
													 print $out,"\n";
												 }
											 }
										 }
								 }
								 }
								else
								{
									print "Annotation Repeat is not present in the entry: $id\n"; 
								}
								close (fw);
unlink ("cont.txt");
}


sub Seqconf {
open (fw,">cont.txt");	
print "Data for Sequence conflict:\n";
$id = shift;
chomp $id;
Geteq($id);
								$feturl="https://www.uniprot.org/uniprot/?query=id:$id&format=tab&columns=feature(SEQUENCE CONFLICT)";
								$mech->get($feturl);
								 $featcontout = $mech->content();
								 @flcont=split("\n",$featcontout);
								foreach $hbk (@flcont) {
									chomp $hbk;
									if ($hbk=~/^\S+/) {
										print fw "$hbk\n";
									}
								}

								$out=`findstr /R /N "^" cont.txt | find /C ":"`;
								if ($out>1) {
								 @actdt=split("\n",$featcontout);
								 foreach $acdtln (@actdt) {
									 chomp $acdtln;									 
									  if ($acdtln=~/^CONFLICT/g) {	
										  $acdtln=~s/\s+CONFLICT/__CONFLICT/g;
											 @splact=split("__",$acdtln);											  
											 foreach $echac (@splact) {
												 $echac=~s/^\s+//g;
												 $echac=~s/\s+$//g;												 
												 if ($echac=~/CONFLICT\s+(\d+)/) {
													 $stactst=$1;
													 $st=$stactst-1;
													 $end=$2;
													 #@dsbd=split("",$seq);
													 $temp = substr($seq, $st,1);
													 $out= ">$id | Sequence conflict | $stactst\n$conf\n";
													 print $out,"\n";
												 }
											 }
										 }
								 }
								 }
								else
								{
									print "Annotation Sequence conflict is not present in the entry: $id\n"; 
								}
								close (fw);
unlink ("cont.txt");
}


sub SeqUnc { 
open (fw,">cont.txt");
print "Data for Sequence uncertainty:\n";
$id = shift;
chomp $id;
Geteq($id);
								$feturl="https://www.uniprot.org/uniprot/?query=id:$id&format=tab&columns=feature(SEQUENCE UNCERTAINTY)";
								$mech->get($feturl);
								 $featcontout = $mech->content();
								 @flcont=split("\n",$featcontout);
								foreach $hbk (@flcont) {
									chomp $hbk;
									if ($hbk=~/^\S+/) {
										print fw "$hbk\n";
									}
								}

								$out=`findstr /R /N "^" cont.txt | find /C ":"`;
								if ($out>1) {
								 @actdt=split("\n",$featcontout);
								 foreach $acdtln (@actdt) {
									 chomp $acdtln;									 
									  if ($acdtln=~/^UNSURE/g) {
										  $acdtln=~s/\s+UNSURE/__UNSURE/g;
											 @splact=split("__",$acdtln);											  
											 foreach $echac (@splact) {
												 $echac=~s/^\s+//g;
												 $echac=~s/\s+$//g;
												 #print $echac,"\n";
												 if ($echac=~/UNSURE\s+(\d+)\;\s+\/note\=\"(.+)\"/) {
													 $stactst=$1;
													 $st=$stactst-1;
													 $end=$1;
													 $unc=$2;
													 $unc=~s/\"\;\s+\/.+//g;
													 $tot=$end-$stactst;
													 $totseq=$tot+1;
													 $temp = substr($seq, $st,$totseq);
													 $out= ">$id | Sequence uncertainty : $unc | $stactst-$end\n$temp\n";
													 print $out,"\n";
												 }
											 }
										 }
								 }
								 }
								else
								{
									print "Annotation Sequence uncertainty is not present in the entry: $id\n"; 
								}
								close (fw);
unlink ("cont.txt");
}


sub Signpep { 
open (fw,">cont.txt");
print "Data for Signal peptide:\n";
$id = shift;
chomp $id;
Geteq($id);
								$feturl="https://www.uniprot.org/uniprot/?query=id:$id&format=tab&columns=feature(SIGNAL)";
								$mech->get($feturl);
								 $featcontout = $mech->content();
								 @flcont=split("\n",$featcontout);
								foreach $hbk (@flcont) {
									chomp $hbk;
									if ($hbk=~/^\S+/) {
										print fw "$hbk\n";
									}
								}

								$out=`findstr /R /N "^" cont.txt | find /C ":"`;
								if ($out>1) {
								 @actdt=split("\n",$featcontout);
								 foreach $acdtln (@actdt) {
									 chomp $acdtln;									 
									  if ($acdtln=~/^SIGNAL/g) {
										  $acdtln=~s/\s+SIGNAL/__SIGNAL/g;
											 @splact=split("__",$acdtln);											  
											 foreach $echac (@splact) {
												 $echac=~s/^\s+//g;
												 $echac=~s/\s+$//g;												 
												 if ($echac=~/SIGNAL\s+(\d+)\.\.(\d+)/) {
													 $stactst=$1;
													 $st=$stactst-1;
													 $end=$2;
													 $tot=$end-$stactst;
													 $totseq=$tot+1;
													 $temp = substr($seq, $st,$totseq);
													 $out= ">$id | Signal peptide | $stactst-$end\n$temp\n";
													 print $out,"\n";
												 }
											 }
										 }
								 }
								 }
								else
								{
									print "Annotation Signal peptide is not present in the entry: $id\n"; 
								}
								close (fw);
unlink ("cont.txt");
}



sub Site { 
open (fw,">cont.txt");
print "Data for Site:\n";
$id = shift;
chomp $id;
Geteq($id);
								$feturl="https://www.uniprot.org/uniprot/?query=id:$id&format=tab&columns=feature(SITE)";
								$mech->get($feturl);
								 $featcontout = $mech->content();
								 @flcont=split("\n",$featcontout);
								foreach $hbk (@flcont) {
									chomp $hbk;
									if ($hbk=~/^\S+/) {
										print fw "$hbk\n";
									}
								}

								$out=`findstr /R /N "^" cont.txt | find /C ":"`;
								if ($out>1) {
								 @actdt=split("\n",$featcontout);
								 foreach $acdtln (@actdt) {
									 chomp $acdtln;									 
									  if ($acdtln=~/^SITE/g) {
										  $acdtln=~s/\s+SITE/__SITE/g;
											 @splact=split("__",$acdtln);											  
											 foreach $echac (@splact) {
												 $echac=~s/^\s+//g;
												 $echac=~s/\s+$//g;												 
												 if ($echac=~/SITE\s+(\d+)\.\.(\d+)\;\s+\/note\=\"(.+)\"/) {
													 $stactst=$1;
													 $st=$stactst-1;
													 $end=$2;
													 $dnm=$3;
													 $dnm=~s/\"\;\s+\/.+//g;
													 $tot=$end-$stactst;
													 $totseq=$tot+1;
													 $temp = substr($seq, $st,$totseq);
													 $out= ">$id | Site: $dnm | $stactst-$end\n$temp\n";
													 print $out,"\n";
												 }
											 }
										 }
								 }
								 }
								else
								{
									print "Annotation Site is not present in the entry: $id\n"; 
								}
								close (fw);
unlink ("cont.txt");
}


sub Topodom { 
open (fw,">cont.txt");
print "Data for Topological domain:\n";
$id = shift;
chomp $id;
Geteq($id);
								$feturl="https://www.uniprot.org/uniprot/?query=id:$id&format=tab&columns=feature(TOPOLOGICAL DOMAIN)";
								$mech->get($feturl);
								 $featcontout = $mech->content();
								 @flcont=split("\n",$featcontout);
								foreach $hbk (@flcont) {
									chomp $hbk;
									if ($hbk=~/^\S+/) {
										print fw "$hbk\n";
									}
								}

								$out=`findstr /R /N "^" cont.txt | find /C ":"`;
								if ($out>1) {
								 @actdt=split("\n",$featcontout);
								 foreach $acdtln (@actdt) {
									 chomp $acdtln;									 
									  if ($acdtln=~/^TOPO_DOM/g) {	
										  $acdtln=~s/\s+TOPO_DOM/__TOPO_DOM/g;
											 @splact=split("__",$acdtln);											  
											 foreach $echac (@splact) {
												 $echac=~s/^\s+//g;
												 $echac=~s/\s+$//g;												 
												 if ($echac=~/TOPO_DOM\s+(\d+)\.\.(\d+)\;\s+\/note\=\"(.+)\"\;/) {
													 $stactst=$1;
													 $st=$stactst-1;
													 $end=$2;
													 $dnm=$3;
													 $dnm=~s/\"\;\s+\/.+//g;
													 $tot=$end-$stactst;
													 $totseq=$tot+1;
													 $temp = substr($seq, $st,$totseq);
													 $out= ">$id | Topological domain: $dnm | $stactst-$end\n$temp\n";
													 print $out,"\n";
												 }
											 }
										 }
								 }
								 }
								else
								{
									print "Annotation Topological Domain is not present in the entry: $id\n"; 
								}
								close (fw);
unlink ("cont.txt");
}


sub Transpep { 
open (fw,">cont.txt");
print "Data for Transit peptide:\n";
$id = shift;
chomp $id;
Geteq($id);
								$feturl="https://www.uniprot.org/uniprot/?query=id:$id&format=tab&columns=feature(TRANSIT)";
								$mech->get($feturl);
								 $featcontout = $mech->content();
								 @flcont=split("\n",$featcontout);
								foreach $hbk (@flcont) {
									chomp $hbk;
									if ($hbk=~/^\S+/) {
										print fw "$hbk\n";
									}
								}

								$out=`findstr /R /N "^" cont.txt | find /C ":"`;
								if ($out>1) {
								 @actdt=split("\n",$featcontout);
								 foreach $acdtln (@actdt) {
									 chomp $acdtln;									 
									  if ($acdtln=~/^TRANSIT/g) {
										  $acdtln=~s/\s+TRANSIT/__TRANSIT/g;
											 @splact=split("__",$acdtln);											  
											 foreach $echac (@splact) {
												 $echac=~s/^\s+//g;
												 $echac=~s/\s+$//g;												 
												 if ($echac=~/TRANSIT\s+(\d+)\.\.(\d+)/) {
													 $stactst=$1;
													 $st=$stactst-1;
													 $end=$2;
													 $tot=$end-$stactst;
													 $totseq=$tot+1;
													 $temp = substr($seq, $st,$totseq);
													 $out= ">$id | Transit peptide | $stactst-$end\n$temp\n";
													 print $out,"\n";
												 }
											 }
										 }
								 }
								}
								else
								{
									print "Annotation Transit peptide is not present in the entry: $id\n"; 
								}
								close (fw);
unlink ("cont.txt");
}


sub Transmembrane { 
open (fw,">cont.txt");
print "Data for Transmembrane:\n";
$id = shift;
chomp $id;
Geteq($id);
								$feturl="https://www.uniprot.org/uniprot/?query=id:$id&format=tab&columns=feature(TRANSMEMBRANE)";
								$mech->get($feturl);
								 $featcontout = $mech->content();
								 @flcont=split("\n",$featcontout);
								foreach $hbk (@flcont) {
									chomp $hbk;
									if ($hbk=~/^\S+/) {
										print fw "$hbk\n";
									}
								}

								$out=`findstr /R /N "^" cont.txt | find /C ":"`;
								if ($out>1) {
								 @actdt=split("\n",$featcontout);
								 foreach $acdtln (@actdt) {
									 chomp $acdtln;									 
									  if ($acdtln=~/^TRANSMEM/g) {
										  $acdtln=~s/\s+TRANSMEM/__TRANSMEM/g;
											 @splact=split("__",$acdtln);											  
											 foreach $echac (@splact) {
												 $echac=~s/^\s+//g;
												 $echac=~s/\s+$//g;												 
												 if ($echac=~/TRANSMEM\s+(\d+)\.\.(\d+)/) {
													 $stactst=$1;
													 $st=$stactst-1;
													 $end=$2;
													 $tot=$end-$stactst;
													 $totseq=$tot+1;
													 $temp = substr($seq, $st,$totseq);
													 $out= ">$id | Transmembrane: $dnm | $stactst-$end\n$temp\n";
													 print $out,"\n";
												 }
											 }
										 }
								 }
								 }
								else
								{
									print "Annotation Transmembrane is not present in the entry: $id\n"; 
								}
								close (fw);
unlink ("cont.txt");
}


sub Turn {
open (fw,">cont.txt");
$id = shift;
chomp $id;
Geteq($id);
print "Data for Turn:\n";
								$feturl="https://www.uniprot.org/uniprot/?query=id:$id&format=tab&columns=feature(TURN)";
								$mech->get($feturl);
								$featcontout = $mech->content();
								@flcont=split("\n",$featcontout);
								foreach $hbk (@flcont) {
									chomp $hbk;
									if ($hbk=~/^\S+/) {
										print fw "$hbk\n";
									}
								}

								$out=`findstr /R /N "^" cont.txt | find /C ":"`;
								if ($out>1) {
								@actdt=split("\n",$featcontout);
								 foreach $acdtln (@actdt) {
									 chomp $acdtln;
									  if ($acdtln=~/^TURN/g) {
										  $acdtln=~s/\s+TURN/__TURN/g;
											 @splact=split("__",$acdtln);
											 foreach $echac (@splact) {
												 $echac=~s/^\s+//g;
												 $echac=~s/\s+$//g;
												 if ($echac=~/TURN\s+(\d+)\.\.(\d+)/g) {													 
													 $stactst=$1;
													 $st=$stactst-1;
													 $end=$2;
													 $tot=$end-$stactst;
													 $totseq=$tot+1;
													 $st=$stactst-1;
													 $temp = substr($seq, $st,$totseq);
													 $out= ">$id | TURN | $stactst-$end\n$temp\n";
													 print $out,"\n";
												 }
											 }
										 }
								 }
								 }
								else
								{
									print "Annotation Turn is not present in the entry: $id\n"; 
								}
								close (fw);
unlink ("cont.txt");
}


sub Zincfing { 
open (fw,">cont.txt");
print "Data for Zinc finger:\n";
$id = shift;
chomp $id;
Geteq($id);
								$feturl="https://www.uniprot.org/uniprot/?query=id:$id&format=tab&columns=feature(ZINC FINGER)";
								$mech->get($feturl);
								 $featcontout = $mech->content();
								 @flcont=split("\n",$featcontout);
								foreach $hbk (@flcont) {
									chomp $hbk;
									if ($hbk=~/^\S+/) {
										print fw "$hbk\n";
									}
								}

								$out=`findstr /R /N "^" cont.txt | find /C ":"`;
								if ($out>1) {
								 @actdt=split("\n",$featcontout);
								 foreach $acdtln (@actdt) {
									 chomp $acdtln;									 
									  if ($acdtln=~/^ZN_FING/g) {
										  $acdtln=~s/\s+ZN_FING/__ZN_FING/g;
											 @splact=split("__",$acdtln);											  
											 foreach $echac (@splact) {
												 $echac=~s/^\s+//g;
												 $echac=~s/\s+$//g;												 
												 if ($echac=~/ZN_FING\s+(\d+)\.\.(\d+)\;\s+\/note\=\"(.+)\"\;/) {
													 $stactst=$1;
													 $st=$stactst-1;
													 $end=$2;
													 $dnm=$3;
													 $dnm=~s/\"\;\s+\/.+//g;
													 $tot=$end-$stactst;
													 $totseq=$tot+1;
													 $temp = substr($seq, $st,$totseq);
													 $out= ">$id | Zinc finger: $dnm | $stactst-$end\n$temp\n";
													 print $out,"\n";
												 }
											 }
										 }
								 }
								 }
								else
								{
									print "Annotation Domain is not present in the entry: $id\n"; 
								}
								close (fw);
unlink ("cont.txt");
}

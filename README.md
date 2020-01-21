# PFSE
CustFeat is Perl module for extracting and downloading sequence/s of protein feature/s (Sequence annotations) in customised manner 

CustFeat is perl module developed at Bioinformatics Centre, Savitribai Phule Pune University (Formerly University of Pune) for extracting the sequence of any user-specified feature/s out of a total of 39 features annotations available in UniprotKB Entry. You can read more details at https://www.uniprot.org/help/sequence_annotation To extract the sequence of feature/s, user just have to provide the feature/s name and UniprotKB ID/s. The output will be sequence of user provided feature/s of corresponding UniprotKB Entry/ies. Here is the feature names and its abbreviation which are need to provide in perl code to download its sequence.

ActSite => Active Site 
AltSeq => Alternative sequence 
BetaStrand => Beta Strand 
BindSite => Binding site 
Chain => Chain 
CoilCoil => Coiled coil 
CompBias => Compositional bias 
Crosslink => Cross-link 
DNABind => DNA binding 
DisulBond => Disulfide bond 
Domain => Domain 
Glyco => Glycosylation 
Helix => Helix 
Initmeth => Initiator methionine 
Intramemb => Intramembrane 
Lipidation => Lipidation 
MetBind => Metal binding 
ModRes => Modified residue 
Motif => Motif 
NatVar => Natural variant 
Nonadres => Non-adjacent residues 
Nonstdres => Non-standard residue 
Nonterres => Non-terminal residue 
Nubind => Nucleotide binding 
Peptide => Peptide 
Propeptide => Propeptide 
Region => Region 
Repeat => Repeat 
SeqUnc => Sequence uncertainty 
Seqconf => Sequence conflict 
Signpep => Signal peptide 
Site => Site 
Topodom => Topological domain 
Transmembrane => Transmembrane 
Transpep => Transit peptide 
Turn => Turn 
Zincfing => Zinc finger

See the sample perl script to see how to provide feature name and UniProt Id/s in perl script. How to install this module?? Instructions:

1 - Download and extract the Zip folder from https://github.com/rajivkarbhal/PFSE
2 - Unzip the folder. Folder contains the following files- CustFeat.pm, Sample1.pl, Sample2.pl and List.txt.
3 - Copy the CustFeat.pm file and keep it in the directory "C:\Perl64\site\lib" (path may vary). Refresh the folder.
4 - *Provide Proxy authentication in CustFeat.pm if your network have proxy connection otherwise skip this step.
5 - Run the Perl script from anywhere on from computer

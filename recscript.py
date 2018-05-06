#!/usr/bin/env python
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation

s = [record for record in SeqIO.parse("start.gb","genbank")][0]
e = [record for record in SeqIO.parse("end.gb","genbank")][0]
# kan = [record for record in SeqIO.parse("Isce_Kan.gb","genbank")][0]

#ISce_Kan cassette hard-coded as the SeqRecord object below
kan = SeqRecord(Seq("TAGGGATAACAGGGTAATCGATTTATTCAACAAAGCCACGTTGTGTCTCAAAATCTCTGATGTTACATTGCACAAGATAAAAATATATCATCATGAACAATAAAACTGTCTGCTTACATAAACAGTAATACAAGGGGTGTTATGAGCCATATTCAACGGGAAACGTCTTGCTCGAGGCCGCGATTAAATTCCAACATGGATGCTGATTTATATGGGTATAAATGGGCTCGCGATAATGTCGGGCAATCAGGTGCGACAATCTATCGATTGTATGGGAAGCCCGATGCGCCAGAGTTGTTTCTGAAACATGGCAAAGGTAGCGTTGCCAATGATGTTACAGATGAGATGGTCAGACTAAACTGGCTGACGGAATTTATGCCTCTTCCGACCATCAAGCATTTTATCCGTACTCCTGATGATGCATGGTTACTCACCACTGCGATCCCCGGGAAAACAGCATTCCAGGTATTAGAAGAATATCCTGATTCAGGTGAAAATATTGTTGATGCGCTGGCAGTGTTCCTGCGCCGGTTGCATTCGATTCCTGTTTGTAATTGTCCTTTTAACAGCGATCGCGTATTTCGTCTCGCTCAGGCGCAATCACGAATGAATAACGGTTTGGTTGATGCGAGTGATTTTGATGACGAGCGTAATGGCTGGCCTGTTGAACAAGTCTGGAAAGAAATGCATAAGCTTTTGCCATTCTCACCGGATTCAGTCGTCACTCATGGTGATTTCTCACTTGATAACCTTATTTTTGACGAGGGGAAATTAATAGGTTGTATTGATGTTGGACGAGTCGGAATCGCAGACCGATACCAGGATCTTGCCATCCTATGGAACTGCCTCGGTGAGTTTTCTCCTTCATTACAGAAACGGCTTTTTCAAAAATATGGTATTGATAATCCTGATATGAATAAATTGCAGTTTCATTTGATGCTCGATGAGTTTTTCTAATCAGAATTGGTTAATTGGTTGTAACACTGGC",IUPAC.ambiguous_dna),id="ISce_Kan",name="ISce_Kan",description="ISce_Kan",
features=[
    SeqFeature(FeatureLocation(start=0,end=988),strand=1,type='insertion_seq',qualifiers={'note':'Kan-Isce insertion'}),    
    SeqFeature(FeatureLocation(start=0,end=24),strand=1,type='primer_bind',qualifiers={'note':'Kan-Isce forward primer'}),
    SeqFeature(FeatureLocation(start=965,end=988),strand=-1,type='primer_bind',qualifiers={'note':'Kan-Isce forward primer'}),
    SeqFeature(FeatureLocation(start=17,end=142),strand=1,type='promoter',qualifiers={'note':'KanaR promoter'}),
    SeqFeature(FeatureLocation(start=141,end=954),strand=1,type='CDS',qualifiers={'note':'aphAI'})
]
)

# print kan.features

# with open('kan.gb','w') as f:
#         SeqIO.write(kan,f,'genbank')

# print len(s)
# print len(e)
# print s.features

def delfind(s,e):
    """Takes a starting genbank seqRecord s and an ending genbank seqRecord e
    with e having a single contiguous deletion relative to s
    Returns the starting and ending string positions (nucleotide position - 1)
    in s that are deleted in e"""

    mutlen = len(s) - len(e)
    pos = 0
    for i in range(len(s)):
        if s.seq[i] != e.seq[i]:
            pos = i
            break

    for x in range(-10,10):
        endseq = s.seq[pos + x + mutlen:]

        if len(s.seq[:pos]) + mutlen + len(endseq) == len(s.seq):
            # print "This value of x = %s fits the length of s" % x
            # print endseq[:10]
            
            if endseq == e.seq[-len(endseq):]:
                # print "The deleted bases in starting seq are bases %s to %s" % (pos+1+x,pos+mutlen+x)
                a = pos + x
                b = pos + mutlen + x
                break
    return a,b

def delFiles(s,e,a,b,kan):
    """Requires four arguments: (1) the starting genbank seqRecord s,
    (2,3) string positions a and b denoting the deletion in s,
    and (4) a genbank seqRecord kan that contains the Isce_kan cassette

    Returns a concatenated genbank seqRecord and file in which the deleted
    sequence is replaced with the kan cassette, flanked by 40 bp direct repeats,
    which are composed of 20bp 5' of the deletion followed by 20bp 3' of the deletion

    Also returns candidate ultramer sequences and the theoretical PCR product
    """

    # print s.seq[:a][-10:]
    # print s.seq[b:b+20]
    # print s.seq[a-20:a]
    # print s.seq[b:][:10]

    INT = (s[:a] + s[b:b+20] + kan + s[a-20:a] + s[b:])
    # print INT
    # print INT.features

    if INT[a+20-40:a+20].seq == INT[a+20+len(kan):a+20+len(kan)+40].seq:
        # print "Direct repeats!"
        dr1start = a+20-40
        dr1end = a+20
        dr2start = a+20+len(kan)
        dr2end = a+20+len(kan)+40

        dr1feat = SeqFeature(FeatureLocation(start=dr1start,end=dr1end),type='misc_feature',qualifiers={'note':"DirectRepeat"})
        dr2feat = SeqFeature(FeatureLocation(start=dr2start,end=dr2end),type='misc_feature',qualifiers={'note':"DirectRepeat"})
        INT.features.append(dr1feat)
        INT.features.append(dr2feat)

    flank1start = a-40
    flank1end = a
    flank2start = a+20+len(kan)+20
    flank2end = a+20+len(kan)+20+40

    homflank1 = SeqFeature(FeatureLocation(start=flank1start,end=flank1end),type='misc_feature',qualifiers={'note':"HomologyFlank1"})
    homflank2 = SeqFeature(FeatureLocation(start=flank2start,end=flank2end),type='misc_feature',qualifiers={'note':"HomologyFlank2"})        

    INT.features.append(homflank1)
    INT.features.append(homflank2)

    FwPrimerEnd = dr1end + 24
    RvPrimerStart = dr2start - 23

    FwPrimer = SeqFeature(FeatureLocation(start=flank1start,end=FwPrimerEnd),type='primer_bind',strand=1,
    qualifiers={'note':"Fw primer"})
    RvPrimer = SeqFeature(FeatureLocation(start=RvPrimerStart,end=flank2end),type='primer_bind',strand=-1,
    qualifiers={'note':"Rv primer"})

    INT.features.append(FwPrimer)
    INT.features.append(RvPrimer)

    with open('INT.gb','w') as f:
        SeqIO.write(INT,f,'genbank')
    
    Fw = INT[flank1start:FwPrimerEnd]
    Rv = INT[RvPrimerStart:flank2end].reverse_complement()

    with open('Fw_primer.gb','w') as f:
        SeqIO.write(Fw,f,'genbank')

    with open('Rv_primer.gb','w') as f:
        SeqIO.write(Rv,f,'genbank')

    PCRprod = INT[flank1start:flank2end]

    with open('PCRprod.gb','w') as f:
        SeqIO.write(PCRprod,f,'genbank')
    
    RES = (INT[:dr1end] + INT[dr2end:])

    if RES.seq == e.seq:
        # print "Success!"
        with open('RES.gb','w') as f:
            SeqIO.write(RES,f,'genbank')


a,b = delfind(s,e)

delFiles(s,e,a,b,kan)
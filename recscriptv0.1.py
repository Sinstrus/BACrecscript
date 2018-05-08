#!/usr/bin/env python
import os
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation
from math import ceil

path = os.path.dirname(os.path.abspath( __file__ ))

with open(path + '/start.gb','r') as f:
    s = [record for record in SeqIO.parse(f,"genbank")][0]

with open(path + '/end.gb','r') as f:
    e = [record for record in SeqIO.parse(f,"genbank")][0]

limit = 90

#ISce_Kan cassette and pSP72 hard-coded as the SeqRecord object below
kan = SeqRecord(Seq("TAGGGATAACAGGGTAATCGATTTATTCAACAAAGCCACGTTGTGTCTCAAAATCTCTGATGTTACATTGCACAAGATAAAAATATATCATCATGAACAATAAAACTGTCTGCTTACATAAACAGTAATACAAGGGGTGTTATGAGCCATATTCAACGGGAAACGTCTTGCTCGAGGCCGCGATTAAATTCCAACATGGATGCTGATTTATATGGGTATAAATGGGCTCGCGATAATGTCGGGCAATCAGGTGCGACAATCTATCGATTGTATGGGAAGCCCGATGCGCCAGAGTTGTTTCTGAAACATGGCAAAGGTAGCGTTGCCAATGATGTTACAGATGAGATGGTCAGACTAAACTGGCTGACGGAATTTATGCCTCTTCCGACCATCAAGCATTTTATCCGTACTCCTGATGATGCATGGTTACTCACCACTGCGATCCCCGGGAAAACAGCATTCCAGGTATTAGAAGAATATCCTGATTCAGGTGAAAATATTGTTGATGCGCTGGCAGTGTTCCTGCGCCGGTTGCATTCGATTCCTGTTTGTAATTGTCCTTTTAACAGCGATCGCGTATTTCGTCTCGCTCAGGCGCAATCACGAATGAATAACGGTTTGGTTGATGCGAGTGATTTTGATGACGAGCGTAATGGCTGGCCTGTTGAACAAGTCTGGAAAGAAATGCATAAGCTTTTGCCATTCTCACCGGATTCAGTCGTCACTCATGGTGATTTCTCACTTGATAACCTTATTTTTGACGAGGGGAAATTAATAGGTTGTATTGATGTTGGACGAGTCGGAATCGCAGACCGATACCAGGATCTTGCCATCCTATGGAACTGCCTCGGTGAGTTTTCTCCTTCATTACAGAAACGGCTTTTTCAAAAATATGGTATTGATAATCCTGATATGAATAAATTGCAGTTTCATTTGATGCTCGATGAGTTTTTCTAATCAGAATTGGTTAATTGGTTGTAACACTGGC",IUPAC.ambiguous_dna),id="ISce_Kan",name="ISce_Kan",description="ISce_Kan",
features=[
    SeqFeature(FeatureLocation(start=0,end=988),strand=1,type='insertion_seq',qualifiers={'note':'Kan-Isce insertion'}),    
    SeqFeature(FeatureLocation(start=0,end=24),strand=1,type='primer_bind',qualifiers={'note':'Kan-Isce forward primer'}),
    SeqFeature(FeatureLocation(start=965,end=988),strand=-1,type='primer_bind',qualifiers={'note':'Kan-Isce reverse primer'}),
    SeqFeature(FeatureLocation(start=17,end=142),strand=1,type='promoter',qualifiers={'note':'KanaR promoter'}),
    SeqFeature(FeatureLocation(start=141,end=954),strand=1,type='CDS',qualifiers={'note':'aphAI'})
]
)

pSP72 = SeqRecord(Seq("GAACTCGAGCAGCTGAAGCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCGGGTACCGAGCTCGAATTCATCGATGATATCAGATCTGCCGGTCTCCCTATAGTGAGTCGTATTAATTTCGATAAGCCAGGTTAACCTGCATTAATGAATCGGCCAACGCGCGGGGAGAGGCGGTTTGCGTATTGGGCGCTCTTCCGCTTCCTCGCTCACTGACTCGCTGCGCTCGGTCGTTCGGCTGCGGCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACAGAATCAGGGGATAACGCAGGAAAGAACATGTGAGCAAAAGGCCAGCAAAAGGCCAGGAACCGTAAAAAGGCCGCGTTGCTGGCGTTTTTCCATAGGCTCCGCCCCCCTGACGAGCATCACAAAAATCGACGCTCAAGTCAGAGGTGGCGAAACCCGACAGGACTATAAAGATACCAGGCGTTTCCCCCTGGAAGCTCCCTCGTGCGCTCTCCTGTTCCGACCCTGCCGCTTACCGGATACCTGTCCGCCTTTCTCCCTTCGGGAAGCGTGGCGCTTTCTCAATGCTCACGCTGTAGGTATCTCAGTTCGGTGTAGGTCGTTCGCTCCAAGCTGGGCTGTGTGCACGAACCCCCCGTTCAGCCCGACCGCTGCGCCTTATCCGGTAACTATCGTCTTGAGTCCAACCCGGTAAGACACGACTTATCGCCACTGGCAGCAGCCACTGGTAACAGGATTAGCAGAGCGAGGTATGTAGGCGGTGCTACAGAGTTCTTGAAGTGGTGGCCTAACTACGGCTACACTAGAAGGACAGTATTTGGTATCTGCGCTCTGCTGAAGCCAGTTACCTTCGGAAAAAGAGTTGGTAGCTCTTGATCCGGCAAACAAACCACCGCTGGTAGCGGTGGTTTTTTTGTTTGCAAGCAGCAGATTACGCGCAGAAAAAAAGGATCTCAAGAAGATCCTTTGATCTTTTCTACGGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGATTTTGGTCATGAGATTATCAAAAAGGATCTTCACCTAGATCCTTTTAAATTAAAAATGAAGTTTTAAATCAATCTAAAGTATATATGAGTAAACTTGGTCTGACAGTTACCAATGCTTAATCAGTGAGGCACCTATCTCAGCGATCTGTCTATTTCGTTCATCCATAGTTGCCTGACTCCCCGTCGTGTAGATAACTACGATACGGGAGGGCTTACCATCTGGCCCCAGTGCTGCAATGATACCGCGAGACCCACGCTCACCGGCTCCAGATTTATCAGCAATAAACCAGCCAGCCGGAAGGGCCGAGCGCAGAAGTGGTCCTGCAACTTTATCCGCCTCCATCCAGTCTATTAATTGTTGCCGGGAAGCTAGAGTAAGTAGTTCGCCAGTTAATAGTTTGCGCAACGTTGTTGCCATTGCTACAGGCATCGTGGTGTCACGCTCGTCGTTTGGTATGGCTTCATTCAGCTCCGGTTCCCAACGATCAAGGCGAGTTACATGATCCCCCATGTTGTGCAAAAAAGCGGTTAGCTCCTTCGGTCCTCCGATCGTTGTCAGAAGTAAGTTGGCCGCAGTGTTATCACTCATGGTTATGGCAGCACTGCATAATTCTCTTACTGTCATGCCATCCGTAAGATGCTTTTCTGTGACTGGTGAGTACTCAACCAAGTCATTCTGAGAATAGTGTATGCGGCGACCGAGTTGCTCTTGCCCGGCGTCAATACGGGATAATACCGCGCCACATAGCAGAACTTTAAAAGTGCTCATCATTGGAAAACGTTCTTCGGGGCGAAAACTCTCAAGGATCTTACCGCTGTTGAGATCCAGTTCGATGTAACCCACTCGTGCACCCAACTGATCTTCAGCATCTTTTACTTTCACCAGCGTTTCTGGGTGAGCAAAAACAGGAAGGCAAAATGCCGCAAAAAAGGGAATAAGGGCGACACGGAAATGTTGAATACTCATACTCTTCCTTTTTCAATATTATTGAAGCATTTATCAGGGTTATTGTCTCATGAGCGGATACATATTTGAATGTATTTAGAAAAATAAACAAATAGGGGTTCCGCGCACATTTCCCCGAAAAGTGCCACCTGACGTCTAAGAAACCATTATTATCATGACATTAACCTATAAAAATAGGCGTATCACGAGGCCCTTTCGTCTCGCGCGTTTCGGTGATGACGGTGAAAACCTCTGACACATGCAGCTCCCGGAGACGGTCACAGCTTGTCTGTAAGCGGATGCCGGGAGCAGACAAGCCCGTCAGGGCGCGTCAGCGGGTGTTGGCGGGTGTCGGGGCTGGCTTAACTATGCGGCATCAGAGCAGATTGTACTGAGAGTGCACCATATGGACATATTGTCGTTAGAACGCGGCTACAATTAATACATAACCTTATGTATCATACACATACGATTTAGGTGACACTATA",IUPAC.ambiguous_dna),id="pSP72",name="pSP72",description="pSP72",
features=[
    SeqFeature(FeatureLocation(start=360,end=980),strand=-1,type='rep_origin',qualifiers={'label':'pBR322_origin'}),    
    SeqFeature(FeatureLocation(start=99,end=118),strand=-1,type='promoter',qualifiers={'label':'T7_promoter'}),
    SeqFeature(FeatureLocation(start=2036,end=2065),strand=-1,type='promoter',qualifiers={'label':'AmpR_promoter'}),
    SeqFeature(FeatureLocation(start=1134,end=1995),strand=-1,type='CDS',qualifiers={'label':'Ampicillin'})
]
)

def output_file(seq_record,filepath):
    """Takes a SeqRecord object and a filepath (assuming the directory exists).
    Outputs a genbank file containing the SeqRecord to the filepath location"""
    with open(filepath,'w') as f:
        SeqIO.write(seq_record,f,'genbank')
    
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
    """Requires five arguments: (1,2) the starting and ending genbank seqRecords s and e,
    (3,4) string positions a and b denoting the deletion in s,
    and (5) a genbank seqRecord kan that contains the Isce_kan cassette

    Returns a concatenated genbank seqRecord and file in which the deleted
    sequence is replaced with the kan cassette, flanked by 40 bp direct repeats,
    which are composed of 20bp 5' of the deletion followed by 20bp 3' of the deletion

    Also returns candidate ultramer sequences and the theoretical PCR product
    """
    INT = (s[:a] + s[b:b+20] + kan + s[a-20:a] + s[b:])

    if INT[a+20-40:a+20].seq == INT[a+20+len(kan):a+20+len(kan)+40].seq:
        dr1start = a+20-40
        dr1end = a+20
        dr2start = a+20+len(kan)
        dr2end = a+20+len(kan)+40

        dr1feat = SeqFeature(FeatureLocation(start=dr1start,end=dr1end),type='misc_feature',qualifiers={'note':"DirectRepeat"})
        dr2feat = SeqFeature(FeatureLocation(start=dr2start,end=dr2end),type='misc_feature',qualifiers={'note':"DirectRepeat"})
        INT.features.append(dr1feat)
        INT.features.append(dr2feat)
    
    RES = (INT[:dr1end] + INT[dr2end:])

    if RES.seq == e.seq:
        print "Success!"

        flank1start = a-40
        flank1end = a
        flank2start = a+20+len(kan)+20
        flank2end = a+20+len(kan)+20+40

        homflank1 = SeqFeature(FeatureLocation(start=flank1start,end=flank1end),type='misc_feature',qualifiers={'note':"HomologyFlank1"})
        homflank2 = SeqFeature(FeatureLocation(start=flank2start,end=flank2end),type='misc_feature',qualifiers={'note':"HomologyFlank2"})

        FwPrimerEnd = dr1end + 24
        RvPrimerStart = dr2start - 23

        FwPrimer = SeqFeature(FeatureLocation(start=flank1start,end=FwPrimerEnd),type='primer_bind',strand=1,qualifiers={'note':"Fw primer"})
        RvPrimer = SeqFeature(FeatureLocation(start=RvPrimerStart,end=flank2end),type='primer_bind',strand=-1,qualifiers={'note':"Rv primer"})

        INT.features.extend((homflank1,homflank2,FwPrimer,RvPrimer))
        
        Fw = INT[flank1start:FwPrimerEnd]
        Rv = INT[RvPrimerStart:flank2end].reverse_complement()

        PCRprod = INT[flank1start:flank2end]

        RES = (INT[:dr1end] + INT[dr2end:]) #reconstruct RES with new features in INT

        output_file(RES,path + '/RES.gb')
        output_file(INT,path + '/INT.gb')
        output_file(Fw,path + '/Fw_primer.gb')
        output_file(Rv,path + '/Rv_primer.gb')
        output_file(PCRprod,path + '/PCRprod.gb')
    else:
        raise ValueError('ERROR: unable to construct sequences! Aborting...')

def med_insfind(s,e):
    """Takes a starting genbank seqRecord s and an ending genbank seqRecord e
    with e having a single contiguous insertion relative to s.
    Returns the string positions (nucleotide position - 1)
    in s, between which the extra sequence has been inserted,
    as well as the inserted sequence itself as a genbank SeqRecord"""

    insertseq = ""
    mutlen = len(e) - len(s)

    while True:
        try:
            with open(path + "/insert.gb",'r') as f:
                insertseq = [record for record in SeqIO.parse(f,"genbank")][0].seq
            
            inslen = len(insertseq)
            if insertseq in e.seq and inslen == mutlen:
                break
            else:
                raise IOError
        except IOError:
            print "Unable to find valid 'insert' genbank file!"
            while True:
                insertseq = raw_input("Paste the insert sequence: ").upper()
                inslen = len(insertseq)
                if insertseq in e.seq:
                    if inslen == mutlen:
                        break
                    else:
                        print "Input sequence found in 'end' sequence, but is too short! Try again..."
                        continue
                else:
                    print "Input sequence not found in 'end' sequence! Try again..."
                    continue
            break
    
    inscount = e.seq.count(insertseq)

    if inscount == 1:
        loc = e.seq.find(insertseq)
        insert = e[loc:loc + inslen]
    else:
        for i in range(len(s)):
            check1 = e.seq[i:i+inslen] == insertseq #checks if the current window matches insertseq
            check2 = e.seq[:i] == s.seq[:i] #checks if e and s have identical prefixes up to the current window
            check3 = e.seq[i+inslen:] == s.seq[-len(e.seq[i+inslen:]):] #checks if e and s have identical suffixes after the current window

            if check1 and check2 and check3:
                # print "Found location of insert!"
                loc = i
                insert = e[loc:loc+inslen]
                break
    
    a = loc
    b = a
    return a,b,insert

def med_insFiles(s,e,a,b,kan,insert):
    """Requires five arguments: (1,2) the starting and ending genbank seqRecords s and e,
    (3,4) string positions a and b denoting the insertion in s,
    (5) a genbank seqRecord kan that contains the Isce_kan cassette,
    and (6) a genbank seqRecord containing the insertion sequence

    Returns a concatenated genbank seqRecord and file in which a kan cassette, flanked by 40 bp direct repeats plus
    the rest of the inserted sequence, is inserted between base positions a and b in sequence s.

    Also returns candidate ultramer sequences and the theoretical PCR product
    """

    ins_len = len(insert)
    # print ins_len

    if ins_len % 2 == 0:
        k = (ins_len - 40)/2 #k is the number of bases on either side of the direct repeat to complete the inserted sequence
        drstart = k
        drend = k + 40
    else:
        k_left = (ins_len - 40)/2 #when the insert length is uneven, same as above, except the left side is 1 bp less than k on the right
        # k_right = ins_len - 40 - k_left

        drstart = k_left
        drend = k_left + 40
    
    drfeat = SeqFeature(FeatureLocation(start=drstart,end=drend),type='misc_feature',qualifiers={'note':"DirectRepeat"})
    leftfeat = SeqFeature(FeatureLocation(start=0,end=drend),type='misc_feature',qualifiers={'note':"insert_leftfrag"})
    rightfeat = SeqFeature(FeatureLocation(start=drstart,end=ins_len),type='misc_feature',qualifiers={'note':"insert_rightfrag"})
    insert.features.extend((drfeat,leftfeat,rightfeat))

    INT = (s[:a] + insert[:drend] + kan + insert[drstart:] + s[b:])

    RES = INT[:a+drend] + INT[a+drend+len(kan)+40:] #constructs the RES SeqRecord

    if RES.seq == e.seq: #Is the RES identical to the end sequence?
        print "Success!"

        flank1start = a - 40
        flank1end = a
        flank2start = a + len(insert[:drend]) + len(kan) + len(insert[drstart:])
        flank2end = flank2start + 40

        homflank1 = SeqFeature(FeatureLocation(start=flank1start,end=flank1end),type='misc_feature',qualifiers={'note':"HomologyFlank1"})
        homflank2 = SeqFeature(FeatureLocation(start=flank2start,end=flank2end),type='misc_feature',qualifiers={'note':"HomologyFlank2"})

        FwPrimerEnd = a + len(insert[:drend]) + 24
        RvPrimerStart = a + len(insert[:drend]) + len(kan) - 23

        FwPrimer = SeqFeature(FeatureLocation(start=flank1start,end=FwPrimerEnd),type='primer_bind',strand=1,qualifiers={'note':"Fw primer"})
        RvPrimer = SeqFeature(FeatureLocation(start=RvPrimerStart,end=flank2end),type='primer_bind',strand=-1,qualifiers={'note':"Rv primer"})

        INT.features.extend((homflank1,homflank2,FwPrimer,RvPrimer))

        Fw = INT[flank1start:FwPrimerEnd]
        Rv = INT[RvPrimerStart:flank2end].reverse_complement()

        RES = INT[:a+drend] + INT[a+drend+len(kan)+40:] #reconstructs the RES SeqRecord

        PCRprod = INT[flank1start:flank2end]

        output_file(RES,path + '/RES.gb')
        output_file(INT,path + '/INT.gb')
        output_file(Fw,path + '/Fw_primer.gb')
        output_file(Rv,path + '/Rv_primer.gb')
        output_file(PCRprod,path + '/PCRprod.gb')
    else:
        raise ValueError('ERROR: unable to construct sequences! Aborting...')

def sm_insFiles(s,e,a,b,kan,insert):
    """Requires five arguments: (1,2) the starting and ending genbank seqRecords s and e,
    (3,4) string positions a and b denoting the insertion in s,
    (5) a genbank seqRecord kan that contains the Isce_kan cassette,
    and (6) a genbank seqRecord containing the insertion sequence

    Returns a concatenated genbank seqRecord and file in which a kan cassette, flanked by 40 bp direct repeats,
    is inserted between base positions a and b in sequence s. Each direct repeat is composed of the inserted sequence
    plus sequence from left and right of the insertion location to get up to 40 bp total.

    Homology flanks are 40 bp, taking the above direct repeat sequence into account.

    Also returns candidate ultramer sequences and the theoretical PCR product
    """

    ins_len = len(insert)
    # print ins_len

    if ins_len % 2 == 0:
        k = (40 - ins_len)/2 #k is the number of bases on either side of the direct repeat to complete the inserted sequence
        # drstart = k
        # drend = k + 40

        k_left = k
        k_right = k

        DRseq = s[a-k:a] + insert + s[b:b+k]

        # print DRseq.seq
    else:
        k_left = (40-ins_len)/2 #when the insert length is uneven, same as above, except the left side is 1 bp less than k on the right
        k_right = 40 - ins_len - k_left

        DRseq = s[a - k_left:a] + insert + s[b:b + k_right]

        # print DRseq.seq
    
    drfeat = SeqFeature(FeatureLocation(start=0,end=len(DRseq)),type='misc_feature',qualifiers={'note':"DirectRepeat"})
    DRseq.features.append(drfeat)

    INT = (s[:a - k_left] + DRseq + kan + DRseq + s[b + k_right:])

    # output_file(INT,path + '/INT.gb')

    RES = INT[:a - k_left + len(DRseq)] + INT[a - k_left + len(DRseq) + len(kan) + len(DRseq):] #constructs the RES SeqRecord

    if RES.seq == e.seq: #Is the RES identical to the end sequence?
        print "Success!"

        flank1start = a - 40
        flank1end = a
        flank2start = a - k_left + len(DRseq) + len(kan) + len(DRseq) - k_right
        flank2end = flank2start + 40

        homflank1 = SeqFeature(FeatureLocation(start=flank1start,end=flank1end),type='misc_feature',qualifiers={'note':"HomologyFlank1"})
        homflank2 = SeqFeature(FeatureLocation(start=flank2start,end=flank2end),type='misc_feature',qualifiers={'note':"HomologyFlank2"})

        FwPrimerEnd = a - k_left + len(DRseq) + 24
        RvPrimerStart = a - k_left + len(DRseq) + len(kan) - 23

        FwPrimer = SeqFeature(FeatureLocation(start=flank1start,end=FwPrimerEnd),type='primer_bind',strand=1,qualifiers={'note':"Fw primer"})
        RvPrimer = SeqFeature(FeatureLocation(start=RvPrimerStart,end=flank2end),type='primer_bind',strand=-1,qualifiers={'note':"Rv primer"})

        INT.features.extend((homflank1,homflank2,FwPrimer,RvPrimer))

        Fw = INT[flank1start:FwPrimerEnd]
        Rv = INT[RvPrimerStart:flank2end].reverse_complement()

        RES = INT[:a - k_left + len(DRseq)] + INT[a - k_left + len(DRseq) + len(kan) + len(DRseq):] #constructs the RES SeqRecord
    
        PCRprod = INT[flank1start:flank2end]

        output_file(RES,path + '/RES.gb')
        output_file(INT,path + '/INT.gb')
        output_file(Fw,path + '/Fw_primer.gb')
        output_file(Rv,path + '/Rv_primer.gb')
        output_file(PCRprod,path + '/PCRprod.gb')
    else:
        raise ValueError('ERROR: unable to construct sequences! Aborting...')

def lg_insFiles(s,e,a,b,kan,pSP72,insert):
    """
    Requires seven arguments: (1,2) the starting and ending genbank seqRecords s and e,
    (3,4) string positions a and b denoting the insertion in s,
    (5,6) a genbank seqRecord kan that contains the Isce_kan cassette and pSP72,
    and (7) a genbank seqRecord containing the insertion sequence

    ...Returns WHAT?

    Also returns candidate ultramer sequences and the theoretical PCR product

    GENEBLOCK PRICES
    125     - 500bp     $89
    501     - 750bp     $129
    751     - 1000bp    $149
    1001    - 1250bp    $209
    1251    - 1500bp    $249
    1501    - 1750bp    $289
    1751    - 2000bp    $329
    and so on...
    """

    ins_len = len(insert)
    # print ins_len

    if ins_len % 2 == 0:
        k = (ins_len - 40)/2 #k is the number of bases on either side of the direct repeat to complete the inserted sequence
        drstart = k
        drend = k + 40
    else:
        k_left = (ins_len - 40)/2 #when the insert length is uneven, same as above, except the left side is 1 bp less than k on the right
        # k_right = ins_len - 40 - k_left

        drstart = k_left
        drend = k_left + 40
    
    drfeat = SeqFeature(FeatureLocation(start=drstart,end=drend),type='misc_feature',qualifiers={'note':"DirectRepeat"})
    leftfeat = SeqFeature(FeatureLocation(start=0,end=drend),type='misc_feature',qualifiers={'note':"insert_leftfrag"})
    rightfeat = SeqFeature(FeatureLocation(start=drstart,end=ins_len),type='misc_feature',qualifiers={'note':"insert_rightfrag"})
    insert.features.extend((drfeat,leftfeat,rightfeat))

    SP72_INT = pSP72[:81] + s[a-40:a] + insert[:drend] + kan + insert[drstart:] + s[b:b+40] + pSP72[81:]

    corelen1 = len(s[a-40:a] + insert[:drend])
    gibslen1 = 30
    checklen1 = (corelen1 + gibslen1*2) % 250

    if checklen1 <= 20:
        gibslen1 = gibslen1 - int(ceil(checklen1/2.0))

    corelen2 = len(insert[drstart:] + s[b:b+40])
    gibslen2 = 30
    checklen2 = (corelen2 + gibslen2*2) % 250

    if checklen2 <= 20:
        gibslen2 = gibslen2 - int(ceil(checklen2/2.0))
    
    FRAG1feat = SeqFeature(FeatureLocation(start=81-gibslen1,end=81+40+drend+gibslen1),type='misc_feature',qualifiers={'note':"FRAG1"})
    SP72_INT.features.append(FRAG1feat)

    output_file(SP72_INT,path + '/SP72_INT.gb')
    

    # INT = (s[:a] + insert[:drend] + kan + insert[drstart:] + s[b:])

    # RES = INT[:a+drend] + INT[a+drend+len(kan)+40:] #constructs the RES SeqRecord

    # if RES.seq == e.seq: #Is the RES identical to the end sequence?
    #     print "Success!"

    #     flank1start = a - 40
    #     flank1end = a
    #     flank2start = a + len(insert[:drend]) + len(kan) + len(insert[drstart:])
    #     flank2end = flank2start + 40

    #     homflank1 = SeqFeature(FeatureLocation(start=flank1start,end=flank1end),type='misc_feature',qualifiers={'note':"HomologyFlank1"})
    #     homflank2 = SeqFeature(FeatureLocation(start=flank2start,end=flank2end),type='misc_feature',qualifiers={'note':"HomologyFlank2"})

    #     FwPrimerEnd = a + len(insert[:drend]) + 24
    #     RvPrimerStart = a + len(insert[:drend]) + len(kan) - 23

    #     FwPrimer = SeqFeature(FeatureLocation(start=flank1start,end=FwPrimerEnd),type='primer_bind',strand=1,qualifiers={'note':"Fw primer"})
    #     RvPrimer = SeqFeature(FeatureLocation(start=RvPrimerStart,end=flank2end),type='primer_bind',strand=-1,qualifiers={'note':"Rv primer"})

    #     INT.features.extend((homflank1,homflank2,FwPrimer,RvPrimer))

    #     Fw = INT[flank1start:FwPrimerEnd]
    #     Rv = INT[RvPrimerStart:flank2end].reverse_complement()

    #     RES = INT[:a+drend] + INT[a+drend+len(kan)+40:] #reconstructs the RES SeqRecord

    #     PCRprod = INT[flank1start:flank2end]

    #     output_file(RES,path + '/RES.gb')
    #     output_file(INT,path + '/INT.gb')
    #     output_file(Fw,path + '/Fw_primer.gb')
    #     output_file(Rv,path + '/Rv_primer.gb')
    #     output_file(PCRprod,path + '/PCRprod.gb')
    # else:
    #     raise ValueError('ERROR: unable to construct sequences! Aborting...')

if len(e) < len(s):
    print "Type of mutation: deletion..."
    print "Calculating positions of deleted bases..."
    a,b = delfind(s,e)
    print "Bases %s to %s have been deleted in starting sequence..." % (a+1,b)
    print "Building ISce_Kan integrate..."
    delFiles(s,e,a,b,kan)
    print "Outputting INT, RES, and primer files..."

elif len(e) - len(s) >= 40 and len(e) - len(s) <= limit:
    print "Type of mutation: medium insertion (between 40 and %s bp in length)..." % limit
    print "Calculating position of insertion..."
    a,b,insert = med_insfind(s,e)
    print "Insertion %s...%s found between bases %s and %s of starting sequence..." % (insert.seq[:20],insert.seq[-20:],a,b+1)    
    med_insFiles(s,e,a,b,kan,insert)

elif len(e) - len(s) < 40 and len(e) - len(s) > 0:
    print "Type of mutation: small insertion (less than 40 bp in length)..."
    print "Calculating position of insertion..."
    a,b,insert = med_insfind(s,e)
    print "Insertion %s found between bases %s and %s of starting sequence..." % (insert.seq,a,b+1)
    sm_insFiles(s,e,a,b,kan,insert)

elif len(e) - len(s) > limit:
    print "Type of mutation: large insertion (greater than %s bp in length)..." % limit
    print "Calculating position of insertion..."
    a,b,insert = med_insfind(s,e)
    print "Insertion %s...%s found between bases %s and %s of starting sequence..." % (insert.seq[:20],insert.seq[-20:],a,b+1)
    lg_insFiles(s,e,a,b,kan,pSP72,insert)

# print pSP72.seq[:81][-20:]
# print pSP72.seq[81:][:20]
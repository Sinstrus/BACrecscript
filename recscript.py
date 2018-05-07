"""
determine the kind of desired mutation

take start and end gb files s and e
compare len
if len(e.seq) > len(s.seq):
    it's an insertion
elif len(e.seq) < len(s.seq):
    it's a deletion
elif len(e.seq) == len(s.seq):
    it's some form of substitution mutation

DELETION (easy)


"""

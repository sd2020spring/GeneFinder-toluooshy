# -*- coding: utf-8 -*-
"""
GENE FINDER MINI PROJECT ONE

@author: Tolulope Oshinowo

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    return ''.join(random.sample(s, len(s)))


def get_complement(nucleotide):
    #Making sure that complements work properly.
    """
    >>> get_complement("A")
    'T'
    >>> get_complement("C")
    'G'
    """
    if(nucleotide == 'A'):
        return 'T'
    elif(nucleotide == 'T'):
        return 'A'
    elif(nucleotide == 'C'):
        return 'G'
    elif(nucleotide == 'G'):
        return 'C'
    else:
        return 'no nucleotide given'

def get_reverse_complement(dna):
    #Making sure that reverse strand operates properly.
    """
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    length = len(dna)
    output = ''
    while length > 0:
        output = output + get_complement(dna[length-1:length])
        length = length - 1
    return output

def rest_of_ORF(dna):
    #Making sure start codon can be identified and the stop can get cut off in sets of 3.
    """
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    if(dna[0:3] == 'ATG'):
        pos = 0
        length = len(dna)
        output = ''
        while length > 0:
            if(dna[pos:pos+3] == 'TAA'):
                return output
            if(dna[pos:pos+3] == 'TAG'):
                return output
            if(dna[pos:pos+3] == 'TGA'):
                return output
            output = output + dna[pos:pos+3]
            pos = pos + 3
            length = length - 3
        return output
    else:
        return 'no start codon found'

def find_all_ORFs_oneframe(dna):
    #Making sure that the previous function can generate a list of codons can be returned.
    """
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    pos = 0
    full = len(dna)
    length = len(dna)
    string = ''
    output = []
    while length > 0:
        string = rest_of_ORF(dna[pos:full])
        if(dna[pos:pos+3] == 'ATG'):
            output.append(string)
        pos = pos + 3
        string = ''
        length = length - 3
    return output

def find_all_ORFs(dna):
    #Making sure that the previous function can be used to parse through the stand for codons.
    """
    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    pos = 0
    full = len(dna)
    length = len(dna)
    string = ''
    output = []
    while length > 0:
        output += find_all_ORFs_oneframe(dna)
        dna = dna[full-(length-1):full]
        length = length - 3
    return output

def find_all_ORFs_both_strands(dna):
    #Making sure virtually all the prevous functions can come together to return ORFs of both DNA strands.
    """
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    reverse = get_reverse_complement(dna)
    output = []
    output += find_all_ORFs(dna)
    output += find_all_ORFs(reverse)
    i = len(output)-1
    while i > -1:
        if output[i] == 'ATG':
            output.remove(output[i])
        i-=1
    return output

def longest_ORF(dna):
    #Making sure the longest ORF can be selected from a giveln list of ORFs.
    """
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    strands = find_all_ORFs_both_strands(dna)
    length = len(strands)
    x = 0
    i = 0
    n = 0
    strand = ''
    for i in strands:
        y = len(strands[n])
        if y>x:
            strand = strands[n]
        n += 1
    return strand

def longest_ORF_noncoding(dna, num_trials):
    #This method can not really be doctested due to randomization.
    import random
    max = 0
    for i in range(num_trials):
        newdna = ''.join(random.sample(dna, len(dna)))
        if len(longest_ORF(newdna)) > max:
            max = len(longest_ORF(newdna))
    return max

def coding_strand_to_AA(dna):
    #Makes sure that amino acids are being properly transcribed.
    """
    >>> coding_strand_to_AA("ATGCGA")
    'MR'
    >>> coding_strand_to_AA("ATGCCCGCTTT")
    'MPA'
    """
    length = len(dna)
    pos = 0
    aa = ''
    while(length > 2):
        aa += aa_table[dna[pos:pos+3]]
        pos += 3
        length -= 3
    return aa

def gene_finder(dna):
    #Makes sure it all comes together
    from load import load_seq
    dna = load_seq("/home/tolu/GeneFinder/data/X73525.fa")
    threshold = longest_ORF_noncoding(dna, 1500)
    output = []
    protein = ''
    orfs = find_all_ORFs_both_strands(dna)
    i = 0
    for i in range(len(orfs)):
        if(len(orfs[i]) >= threshold):
            protein = coding_strand_to_AA(orfs[i])
            output.append(protein)
    return output

if __name__ == "__main__":
    import doctest
    print(gene_finder(''))
    '''
    doctest.run_docstring_examples(get_complement, globals(), verbose=True)
    doctest.run_docstring_examples(get_reverse_complement, globals(), verbose=True)
    doctest.run_docstring_examples(rest_of_ORF, globals(), verbose=True)
    doctest.run_docstring_examples(find_all_ORFs_oneframe, globals(), verbose=True)
    doctest.run_docstring_examples(find_all_ORFs, globals(), verbose=True)
    doctest.run_docstring_examples(find_all_ORFs_both_strands, globals(), verbose=True)
    doctest.run_docstring_examples(longest_ORF, globals(), verbose=True)
    #doctest.run_docstring_examples(longest_ORF_noncoding, globals(), verbose=True) randomization removes necessity of this test
    doctest.run_docstring_examples(coding_strand_to_AA, globals(), verbose=True)
    '''

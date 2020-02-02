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
    pass


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
    pass


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
            elif(dna[pos:pos+3] == 'TAG'):
                return output
            elif(dna[pos:pos+3] == 'TGA'):
                return output
            output = output + dna[pos:pos+3]
            pos = pos + 3
            length = length - 1
        return output
    else:
        return 'no start codon found'
    pass


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
    pass


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
        length = length - 6
    return output
    pass


def find_all_ORFs_both_strands(dna):
    #Making sure virtually all the prevous functions can come together to return ORFs of both DNA strands.
    """
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    reverse = get_reverse_complement(dna)
    output = []
    marker = 0
    pos = 0
    pos2 = 0
    full = len(dna)
    length = len(dna)
    while length > 0:
        if(dna[pos:pos+3] == 'ATG'):
            if(len(dna[pos:full]) > 6):
                while pos2 < full-3:
                    if(dna[pos2:pos2+3] == 'TAA'):
                        output += find_all_ORFs(dna)[0:(len(find_all_ORFs(dna))-1)]
                    if(dna[pos2:pos2+3] == 'TAG'):
                        output += find_all_ORFs(dna)[0:(len(find_all_ORFs(dna))-1)]
                    if(dna[pos2:pos2+3] == 'TGA'):
                        output += find_all_ORFs(dna)[0:(len(find_all_ORFs(dna))-1)]
                    pos2 = pos2 + 3
            pos = pos + 3
        pos = pos + 1
        length = length - 1
    pos = 0
    pos2 = 0
    full = len(reverse)
    length = len(reverse)
    while length > 0:
        if(reverse[pos:pos+3] == 'ATG'):
            if(len(reverse[pos:full]) > 6):
                while pos2 < full-3:
                    if(reverse[pos2:pos2+3] == 'TAA'):
                        output += find_all_ORFs(reverse)[0:(len(find_all_ORFs(reverse))-1)]
                    if(reverse[pos2:pos2+3] == 'TAG'):
                        output += find_all_ORFs(reverse)[0:(len(find_all_ORFs(reverse))-1)]
                    if(reverse[pos2:pos2+3] == 'TGA'):
                        output += find_all_ORFs(reverse)[0:(len(find_all_ORFs(reverse))-1)]
                    pos2 = pos2 + 3
                output.append(reverse[pos:full])
            else:
                output.append('ATG')
            pos = pos + 3
        pos = pos + 1
        length = length - 1
    return output
    pass


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
        >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    # TODO: implement this
    pass


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    # TODO: implement this
    pass


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    # TODO: implement this
    pass


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # TODO: implement this
    pass

if __name__ == "__main__":
    import doctest
    doctest.run_docstring_examples(get_complement, globals(), verbose=True)
    doctest.run_docstring_examples(get_reverse_complement, globals(), verbose=True)
    doctest.run_docstring_examples(rest_of_ORF, globals(), verbose=True)
    doctest.run_docstring_examples(find_all_ORFs_oneframe, globals(), verbose=True)
    doctest.run_docstring_examples(find_all_ORFs, globals(), verbose=True)
    doctest.run_docstring_examples(find_all_ORFs_both_strands, globals(), verbose=True)

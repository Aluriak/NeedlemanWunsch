# -*- coding: utf-8 -*-
#########################
#       MAIN            
#########################
"""
Proof of concept : Needleman and Wunsch algorithm by recursive dictionnary strategy for sequences alignment.
get methods and example on Bioinformatique et Génomique course (see example function below)
There is many (and better) more other way to do that. (notabily without recursive dictionnary strategy)

usage:
    main.py seqs idt sub alpha beta

with:
    seqs : name of a fasta formatted file where the two first sequences will be used as data
    idt  : score for identity between sequences
    sub  : score for substitution between sequences
    alpha: cost for GAP opening
    beta : cost for GAP extension
"""


#########################
# IMPORTS               #
#########################
import sys
from collections import OrderedDict # for user input





#########################
# NEEDLEMANWUNSCH CLASS #
#########################
class NeedlemanWunsch(dict):
    """
    Dictionnary that implement sequence alignement in the Needleman and Wunsch Algorithm Way.
    There is many way to improve this object !
    Keys are tuple of two integers (i, j) where i describes index in the first sequence
    and j the index in the second one.
    Values are tuple of number and tuple (s, p) where s is the score of (i, j) and p is another key of
    the same object instance, used for get alignment sequences by referencing the previous (i, j).

    __getitem__ is recursive, and compute only what is needed.
    alignment return the pure result of Needleman and Wunsch algorithm.
    """

    # CONSTRUCTOR #################################################################
    def __init__(self, seq1, seq2, score_idty, score_subs, score_gap_open, score_gap_extd):
        """
        @seq1 str first sequence to be compared
        @seq2 str second sequence
        @score_idty number score when seqx match together
        @score_subs number score when seqx don't match
        @score_gap_open number score for gap opening (alpha in equations)
        @score_gap_extd number score for gap extending (beta in equations)
        @return new NeedlemanWunsch object instance, ready to use
        """
        super(dict, self).__init__()
        self.seq1, self.seq2 = seq1, seq2
        self.score_idty, self.score_subs = float(score_idty), float(score_subs)
        self.score_gap_open, self.score_gap_extd = float(score_gap_open), float(score_gap_extd)

        # Initial case of recursion
        self[0,0] = self.identityScore((0,0)), (-1, -1)



    # PUBLIC METHODS ##############################################################
    def identityScore(self, key):
        """
        @param key valable key (i,j) in comparison matrix
        @return score as number for identity or substitution at given coords in sequences
        """
        i, j = key
        return self.score_idty if self.seq1[i] == self.seq2[j] else self.score_subs


    def gapCost(self, gamma):
        """
        @param gamma parameter of function that compute the GAP cost
        @return score number for given GAP length
        """
        return self.score_gap_open + self.score_gap_extd*(gamma-1)



    def alignment(self, key=None):
        """
        @param key the index where alignment research begin. The greater by default.
        @return str list, with the two sequences aligned
        """
        if key is None:
            key = len(self.seq1)-1, len(self.seq2)-1
        aligned1, aligned2 = "", "" # aligned sequences

        # aligned sequences construction
        # walk into the matrix in reverse way, following the reference of each score
        while key in self:
            i, j = key
            n, m = self[key][1] # previous square

            # diagonal: letters are identics
            if i-n == 1 and j-m == 1:
                aligned1 += self.seq1[i]
                aligned2 += self.seq2[j]
            # GAP case
            else:
                if i-n > 0:
                    aligned1 += self.seq1[i:n:-1]
                    aligned2 += "-" * (i-n)
                if j-m > 0:
                    aligned1 += "-" * (j-m)
                    aligned2 += self.seq2[j:m:-1]

            key = n, m


        return [aligned1[::-1], aligned2[::-1]]




    # ACCESSORS ###################################################################
    def __getitem__(self, key):
        """
        @param tuple (i,j) integers index in matrix of comparison
        @return tuple (score, prev), score of square at key coords and way to previous square (cf object doc)
        @note if returned value is None, there is a problem
        """
        i, j = key
        item = None # problem value

        # out of bounds cases
        if i < 0 or j < 0 or i >= len(self.seq1) or j >= len(self.seq2):
            item = None # invalid value
        # already computed cases
        elif key in self:
            item = self.get(key)
            
        # to be computed cases
        else:
            # get max values for diagonal, horizontal and vertical
            # diagonal (only if not on border col or row)
            if i == 0 or j == 0:
                value_score = None 
            else:
                value_score = self[i-1, j-1][0] + self.identityScore(key) 
            self[key] = value_score, (i-1,j-1) # diagonal is the better… at this time
            # horizontal
            for n in range(i):
                col_value = self[n,j][0] - self.gapCost(i-n)
                # if a better score is find
                if value_score is None or col_value > value_score: 
                    value_score = col_value
                    self[key] = value_score, (n,j) # this horizontal is better… at this time
            # vertical 
            for m in range(j):
                row_value = self[i,m][0] - self.gapCost(j-m)
                # if a better score is find
                if value_score is None or row_value > value_score: 
                    value_score = row_value
                    self[key] = value_score, (i,m) # this vertical is better… at this time
            
            # the better score for self[key] is now found
            item = self[key]


        return item


    # CONVERSION ##################################################################
    def __str__(self):
        """
        @return a matrix with scores directly printable
        """
        sself = " |" + "\t|".join([str(_) for _ in self.seq1])
        for j in range(len(self.seq2)):
            sself += "\n" + self.seq2[j] + "|"
            sself += "\t|".join(str(self[i,j][0]) for i in range(len(self.seq1)))
        return sself










#########################
# FUNCTIONS             #
#########################
def example():
    """
    Basic example of use of NeedlemanWunsch object with values of slide 49 of Bioinformatique et Génomique course
    """
    d = NeedlemanWunsch("TACAGTATA", "TATATA", 5, -4, 10, 0.5)

    print("d[1, 1] = " + str(d[1, 1]))
    print("d[3, 1] = " + str(d[3, 1]))
    print("d[4, 1] = " + str(d[4, 1]))
    print("d[7, 3] = " + str(d[7, 3]))
    print("d[8, 5] = " + str(d[8, 5]))
    print(d)
    print("\n\nAlignment sequences :")
    print("\n".join(d.alignment()))



    # unit tests
    assert(d[3, 1][0] == -0.5) 
    assert(d[4, 1][0] == -1.0)
    assert(d[7, 3][0] == -0.5)
    assert(d[8, 5][0] == 19.0)
    assert(d.alignment() == ['TACAGTATA', 'TA---TATA'])
    print("\nUnit tests said that all is ok !")



def readFASTA(text, results=dict()):
    """
    @param text in FASTA format
    @param results dict where results will be put. A new dict by default.
    @return dictionnary with new entry like {sequence_name: sequence_str}
    @note call this function only if biopython can't be used instead
    """
    string = ''
    name = ''
    for line in text.split('\n'):
        if len(line) > 0:
            if line[0] == '>':
                if len(string) > 0:
                    # append a copy of dna in dnas list
                    results[name] = string
                    string = ''
                name = line[1:]
            elif line[0] != '>': 
                # add line without \n last char to current readed dna
                string += line
    # add last line encountered
    if len(string) > 0:
        # append a copy of dna in dnas list
        results[name] = string
    # end
    return results








#########################
# MAIN                  #
#########################
if __name__ == '__main__':
    if len(sys.argv) < 6:
        print(__doc__)
        print("\n########################\n\n" + example.__doc__)
        example() # run something 
    else:
        with open(sys.argv[1], "r") as f:
            data_fasta = f.read()
            seqs = []
            for seq in readFASTA(data_fasta, OrderedDict()).values():
                seqs.append(seq)
            nw = NeedlemanWunsch(seqs[0], seqs[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
            print(nw)
            print("\n".join(nw.alignment()))







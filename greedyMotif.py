

"""
Definitions of subroutines (and subroutines of subroutines ^_^) and their outcomes:

Count(Motifs): Given a collection of Motifs, count the occurrences of each nucleotide in 
Motifs (in this case a collection of 3-mers taken) at each position i in a string of Motifs 
of length j. Return this count as a dictionary.

Profile(Motifs): Given a collection of Motifs, compute the count, and then find the fractional 
probability of encountering a given nucleotide in Motifs at each position i in a string of Motifs of length j. 
Return this profile as a dictionary.

Consensus(Motifs): Given a collection of Motifs, find the most probable sequence of nucleotides in Motifs for 
each given position. Return this consensus sequence as a string.

Score(Motifs): Compute the consensus string, and then summation count all mismatches of this consensus string 
to the strings in Motifs. Return this count as an integer. 

Pr(Text,Profile): Given string 'Text' (in this case one string in Dna) and the probabilities in Profile 
(see above), find the total probability (independent events) of encountering that string. Return this probability as a float.

ProfileMostProbablePattern(Text,Profile,k) or ProfileMostProbablekmer(Text,Profile,k) : Given string 'Text', probabilities 
in Profile and the length of the k-mer, generate the total probabilities of finding a kmer of length k in 'Text' and find 
the k-mer with the highest probability. Return this k-mer as a string. 
"""

text = 'TTACCATGGGACCGCTGACTGATTTCTGGCGTCAGCGTGATGCTGGTGTGGATGACATTCCGGTGCGCTTTGTAAGCAGAGTTTA'
k = 12
profile = {'A': [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.1, 0.2, 0.3, 0.4, 0.5], 
'C': [0.3, 0.2, 0.1, 0.1, 0.2, 0.1, 0.1, 0.4, 0.3, 0.2, 0.2, 0.1], 
'G': [0.2, 0.1, 0.4, 0.3, 0.1, 0.1, 0.1, 0.3, 0.1, 0.1, 0.2, 0.1], 
'T': [0.3, 0.4, 0.1, 0.1, 0.1, 0.1, 0.0, 0.2, 0.4, 0.4, 0.2, 0.3]}



def Count(Motifs):
    k = len(Motifs[0])
    nuc_freq = {symbol: [0]*k for symbol in "ACGT"}
    

    for row in Motifs:
        for index, char in enumerate(row):
            nuc_freq[char][index] += 1
    
    return nuc_freq

def Profile(Motifs):
    k = len(Motifs[0])
    nuc_freq = {symbol: [0]*k for symbol in "ACGT"}
    for row in Motifs:
        for index, char in enumerate(row):
            nuc_freq[char][index] += 1/k    
    return nuc_freq

def Consensus(Motifs):
    k = len(Motifs[0])
    count = Count(Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus

def Score(Motifs):
    consensus = Consensus(Motifs)
    count = 0
    for motif in Motifs:
        for index in range(len(motif)):
            if motif[index] != consensus[index]:
                count += 1
    return count

#This helps to find the probability of a candidate motif given the motif sequence and the matrix profile
def Pr(Text, Profile):
    length = len(Text)
    p = 1
    for i in range(length):
        p = p * Profile[Text[i]][i]        
    return p

def ProfileMostProbableKmer(text, k, profile):      
    bestMer = {"probability":Pr(text[:k], profile), "mer":text[:k]}

    for i in range(1, len(text) - k + 1):           #Range starts from second position of text to the end minus k 
        
        mer = text[i:i+k]                       #Get mer
        pr = Pr(mer, profile)                   #Get probability of mer 
        if  pr > bestMer["probability"]:        #Compare probability of mer to best mer
            bestMer["probability"] = pr
            bestMer["mer"] = mer  

    return bestMer["mer"]       


print(ProfileMostProbableKmer(text, k, profile))

# Input:  A list of kmers Dna, and integers k and t (where t is the number of kmers in Dna)
# Output: GreedyMotifSearch(Dna, k, t)
def GreedyMotifSearch(Dna, k, t):
    BestMotifs = []
    for i in range(t): #range over strings in Dna
        BestMotifs.append(Dna[i][0:k]) # sets best motifs to first 3-mer in each dna string
    n = len(Dna[0])
    for i in range(n-k+1): #range over length of a string in Dna (all of them are the same length)
        Motifs = [] #empty list of Motifs
        Motifs.append(Dna[0][i:i+k]) #get Motif from first string in Dna consecutive possible 3-mers
        for j in range(1,t): #range over from Dna[2nd string] until the end of Dna strings
            P = Profile(Motifs[0:j]) #Get profile of probabilities for Motif previously obtained against the strings below  
            Motifs.append(ProfileMostProbableKmer(Dna[j],k,P)) #append most probable k-mers found in Dna strings 2 to 5 to the list based on the motif obtained from Dna string 1
        if Score(Motifs) < Score(BestMotifs): # create consensus string, and then compute dissimilarity scores for current Motifs and previous BestMotifs, if the current is less (i.e. better) than the previous. 
            BestMotifs = Motifs #update the best as the current Motifs.
    return BestMotifs
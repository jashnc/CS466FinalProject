# Implement Step 2, the algorithm

# Input:
	# (a) A file containing the sequences called "sequences.fa"
	# (b) A file containing the motif length called "motiflength.txt"

# Output
	# (a) A file containing the predicted locations of sites called "predictedsites.txt"
	# (b) A file containing the predicted motif in a specified format called "predictedmotif.txt"

import math

####################################################################################################
#	        										   #
# 					Functions						   #
#												   #
####################################################################################################

####################################################################################################
#       informationContent
#
#       Computes the info content score of an alignment
#       Input: A list of strings (each string is a nucleotide subsequence)
#
#       Output: The information Content score of this alignment
def informationContent(theSequences):
        total = len(theSequences)
        q_B = 1 / total
        
        #make a 2D list (matrix) to represent the profile matrix profileMatrix[col][row]
        #initialize each entry to 0
        profileMatrix = []
        for i in range(len(theSequences[0])):
                newDim = []
                for j in range(4):
                        newDim.append(0)
                profileMatrix.append(newDim)

        # Loop over every sequence provided to populate the profileMatrix
        # increment the proper entries as you loop over them
        for oneSeq in theSequences:
                #print(oneSeq)
                for i in range(len(oneSeq)):
                        #print(i)
                        if(oneSeq[i]=='A'):
                                profileMatrix[i][0] = profileMatrix[i][0] + 1
                        elif(oneSeq[i]=='C'):
                                profileMatrix[i][1] = profileMatrix[i][1] + 1
                        elif(oneSeq[i]=='G'):
                                profileMatrix[i][2] = profileMatrix[i][2] + 1
                        elif(oneSeq[i]=='T'):
                                profileMatrix[i][3] = profileMatrix[i][3] + 1

        #print(profileMatrix)
        #now actually calculate the information content of the profile matrix
        infoContent = 0.0
        for col in profileMatrix:
                for base in col:
                        if(base > 0):     
                                W_BK = base / total
                                newInfo = ( W_BK * math.log( base, 2) )
                                infoContent = infoContent + newInfo

        #return the score!
        return infoContent
#
####################################################################################################

####################################################################################################
#									        		   #
# 			                0) Read in the input files	                           #
#												   #
####################################################################################################

# Read in the motif length
Motif = open('fake_data/motiflength.txt')
motifLength = int(Motif.read())
Motif.close()


# Read in the sequences
Seq = open('fake_data/sequences.fa', 'r')
idx = -1
sequence = []   #list of sequences 
for line in Seq:
        line.rstrip('\r\n')
        if(line[0] == ">"):
                #print("NEW SEQUENCE")
                #print(" ")
                nextSeq = []
                sequence.append(nextSeq)
                idx = idx + 1
        else:
                sequence[idx].append(line.rstrip('\r\n'))
                #print (" ")
                #print (line.rstrip('\r\n'))
                #print (" ")

print(sequence)
Seq.close()



########################################################################################################
#																									   #
# 		                 1) Align sequences 1 and 2 together by trying 								   
#				           every combination and maximizing the score	   							   
#																									   #
########################################################################################################

# Loop over every possible starting positions of seq 1 and 2
# Keep Track of the best alignment and the score for that alignment
bestI = 0
bestJ = 0
bestScore = 0
for i in range(len(sequence[0][0]) - motifLength):
        for j in range(len(sequence[0][0]) - motifLength):
                infoContentInput = []
                infoContentInput.append(sequence[0][0][i:(i+motifLength)])
                infoContentInput.append(sequence[1][0][j:(j+motifLength)])
                if(informationContent(infoContentInput) > bestScore):
                        bestScore = informationContent(infoContentInput)
                        bestI = i
                        bestJ = j

print("Best Alignment of first 2 sequences = ")
print(sequence[0][0][bestI:(bestI+motifLength)])
print(sequence[1][0][bestJ:(bestJ+motifLength)])


########################################################################################################
#																									   #
# 		              1) Iteratively add another sequence to the current alignment by 				   
#				     trying each position in the new sequence and maximizing the score   			   
#																									   #
########################################################################################################






# Implement Step 2, the algorithm

# Input:
        # (a) A file containing the sequences called "sequences.fa"
        # (b) A file containing the motif length called "motiflength.txt"

# Output:
        # (a) A file containing the predicted locations of sites called "predictedsites.txt"
        # (b) A file containing the predicted motif in a specified format called "predictedmotif.txt"

import math
import copy
import sys
import timeit

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
                for i in range(len(oneSeq)):
                        if(oneSeq[i]=='A'):
                                profileMatrix[i][0] = profileMatrix[i][0] + 1
                        elif(oneSeq[i]=='C'):
                                profileMatrix[i][1] = profileMatrix[i][1] + 1
                        elif(oneSeq[i]=='G'):
                                profileMatrix[i][2] = profileMatrix[i][2] + 1
                        elif(oneSeq[i]=='T'):
                                profileMatrix[i][3] = profileMatrix[i][3] + 1

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

def findMotif(Dir):
        # Read in the command line input
        #print(sys.argv[1])
        #Dir = sys.argv[1];

        # Read in the motif length
        Motif = open(Dir + '/motiflength.txt')
        motifLength = int(Motif.read())
        Motif.close()


        # Read in the sequences
        Seq = open(Dir + '/sequences.fa', 'r')
        idx = -1
        sequence = []   #list of sequences 
        for line in Seq:
                line.rstrip('\r\n')
                if(line[0] == ">"):
                        nextSeq = []
                        sequence.append(nextSeq)
                        idx = idx + 1
                else:
                        sequence[idx].append(line.rstrip('\r\n'))

        Seq.close()

        ####################################################################################################
        #                                                                                                  #
        #                                2) Align sequences 1 and 2 together                               #
        #                                                                                                  #
        ####################################################################################################
        # Try every combination and maximize the score to obtain this alignment

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
                        score = informationContent(infoContentInput) 
                        if(score > bestScore):
                                bestScore = score
                                bestI = i
                                bestJ = j

        #keep a variable for the best alignment, add on to it later incrementally
        alignment = []
        alignment.append(sequence[0][0][bestI:(bestI+motifLength)])
        alignment.append(sequence[1][0][bestJ:(bestJ+motifLength)])

        #keep variable for predictedSites
        predictedSites = []
        predictedSites.append(bestI)
        predictedSites.append(bestJ)

        ####################################################################################################
        #                                                                                                  #
        #                              3) Iteratively add the remaining sequences                          #
        #                                                                                                  #
        ####################################################################################################
        # Maintain the current alignment, and try each position of the next sequence to maximize the score 
        #       and align the sequence

        #iterate over all remaining sequences
        for i in range(2,len(sequence)):
                #initialize maximizing variables
                bestJ = 0
                bestScore = 0
                #for this sequence iterate over every possible starting position
                for j in range(len(sequence[0][0]) - motifLength):
                        newAlignment = copy.deepcopy(alignment)
                        newAlignment.append(sequence[i][0][j:(j+motifLength)])
                        score = informationContent(newAlignment)
                        if(score > bestScore):
                                bestScore = score
                                bestJ = j
                #after done iterating, update current alignment and sites
                alignment.append(sequence[i][0][bestJ:(bestJ+motifLength)])
                predictedSites.append(bestJ)

        ####################################################################################################
        #                                                                                                  #
        #                        3) Output the motif we found                                              #
        #                                                                                                  #
        ####################################################################################################

        # Convert the alignment to a profilematrix
        #make a 2D list (matrix) to represent the profile matrix profileMatrix[col][row]
        #initialize each entry to 0
        profileMatrix = []
        for i in range(len(alignment[0])):
                newDim = []
                for j in range(4):
                        newDim.append(0)
                profileMatrix.append(newDim)

        # Loop over every sequence provided to populate the profileMatrix
        # increment the proper entries as you loop over them
        for oneSeq in alignment:
                for i in range(len(oneSeq)):
                        if(oneSeq[i]=='A'):
                                profileMatrix[i][0] = profileMatrix[i][0] + 1
                        elif(oneSeq[i]=='C'):
                                profileMatrix[i][1] = profileMatrix[i][1] + 1
                        elif(oneSeq[i]=='G'):
                                profileMatrix[i][2] = profileMatrix[i][2] + 1
                        elif(oneSeq[i]=='T'):
                                profileMatrix[i][3] = profileMatrix[i][3] + 1

        # Finally, output this profile Matrix
        f = open(Dir + "/predictedmotif.txt", "w")
        f.write(">MOTIF1 " + str(motifLength) + "\n")
        for col in profileMatrix:
                f.write(str(col[0]) + " " + str(col[1]) + " " + str(col[2]) + " " + str(col[3]) + "\n")
        f.write("<\n")
        f.close()

                
        ####################################################################################################
        #                                                                                                  #
        #                        4) Output the predictedSites we found                                     #
        #                                                                                                  #
        ####################################################################################################
        f = open(Dir + "/predictedsites.txt", "w")
        for site in predictedSites:
                f.write("Start Index: " + str(site+1) + " End Index: " + str(site + motifLength) + "\n")

        f.close()


for i in range(1,7+1):       
        for j in range(1,10+1):
                start_time = timeit.default_timer()
                findMotif("data_set/set" + str(i) + "/data" + str(j))
                print("runtime for dataset %d, data %d: %.15f" % (i, j, timeit.default_timer() - start_time))




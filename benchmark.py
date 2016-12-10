# Implement Step 1, generate the datasets

# Input:
	# (a) A positive number called ICPC (“information content per column”)
	# (b) A positive integer called ML (“motif length”)
	# (c) A positive integer called SL (“sequence length”)
	# (d) A positive integer called SC (“sequence count”)

# Output
	# (a) A file containing the sequences called "sequences.fa"
	# (b) A file containing the locations of planted sites called "sites.txt"
	# (c) A file containing the motif in a specified format called "motif.txt"
	# (d) A file containing the motif length called "motiflength.txt"

import random, math, copy

ICPC = 2
ML = 8
SL = 500
SC = 10

#Generate random nucleotide sequences according to input parameters
def generate_sequence(SC, SL):
    nucleotides = ['A','C','G','T']
    sequences = []
    for i in range(0, SC):
        sequence = []
        for j in range(0, SL):
            sequence.append(random.choice(nucleotides))

        sequences.append(sequence)

    return sequences



#Generate random motif of length ML
def generate_motif(SC, SL, ML):
    q_beta = .25

    sample_sequence = generate_sequence(SC, ML)
    nucleotides = ['A','C','G','T']
    profile_matrix = []



    for h in range(ML):
        col = [row[h] for row in sample_sequence]
        row_count = count(col, SC, ICPC)
        profile_matrix.append(row_count)

    pwm = list(zip(*profile_matrix[::-1]))



    return pwm

def count(sequence, SC, ICPC):
    counter = [0,0,0,0]
    for i in range(len(sequence)):
        if sequence[i] == 'A':
            counter[0]+=1
        elif sequence[i] == 'C':
            counter[1]+=1
        elif sequence[i] == 'G':
            counter[2]+=1
        elif sequence[i] == 'T':
            counter[3]+=1

    counter = [(x / SC) for x in counter]
    if counter[0] != 0.0:
        counter[0] = -math.log(counter[0], 2)
    if counter[1] != 0.0:
        counter[1] = -math.log(counter[1], 2)
    if counter[2] != 0.0:
        counter[2] = -math.log(counter[2], 2)
    if counter[3] != 0.0:
        counter[3] = -math.log(counter[3], 2)

    total = ICPC/sum(counter)

    counter = [x * total for x in counter]

    return counter

def generate_binding_sites(pwm, SC):
    binding_sites = []

    for i in range(SC):
        binding_sites.append(generate_binding_site(pwm, SC, SL, ML))

    return binding_sites

def generate_binding_site(pwm, SC, SL, ML):
    bases = ['A', 'C', 'G', 'T']
    binding_site = ''

    print(pwm)
    for i in range(ML):
        random_y = random.randint(0, 3)
        random_x = random.randint(0, ML-1)

        num = math.floor((pwm[random_y][random_x] * random.randint(0, 4))) % 4
        if(random.random() < .5):
            binding_site += bases[num]



    return binding_site

def plant(sequences, binding_sites, SL):

    curr = 0

    temp = []
    for row in sequences:
        rand_starting_index = random.randint(0, SL-ML)
        collapsed_row = ''.join(row)
        print(rand_starting_index)
        print(binding_sites[curr])
        temp2 = copy.deepcopy(collapsed_row)

        #collapsed_row.replace(collapsed_row[rand_starting_index:(rand_starting_index+ML)], binding_sites[curr])
        collapsed_row = collapsed_row[:rand_starting_index] + binding_sites[curr] + collapsed_row[rand_starting_index+len(binding_sites[curr]):]
        print(temp2 == collapsed_row)
        curr += 1
        temp.append(collapsed_row)



    new_sequences = []
    for row in temp:

        new_sequences.append(list(row))


    print(len(new_sequences))
    return new_sequences



sequences = generate_sequence(SC, SL)
print(sequences)
pwm = generate_motif(SC, SL, ML)
binding_sites = generate_binding_sites(pwm, SC)
planted = plant(sequences, binding_sites, SL)

print(sequences)
print(planted)
print(sequences == planted)
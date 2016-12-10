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

import random, math, copy, os


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
def generate_motif(SC, SL, ML, ICPC):
    q_beta = .25

    sample_sequence = generate_sequence(SC, ML)
    nucleotides = ['A','C','G','T']
    profile_matrix = []
    for h in range(ML):
        col = [row[h] for row in sample_sequence]
        row_count = count(col, SC, ICPC, ML)
        profile_matrix.append(row_count)
    pwm = list(zip(*profile_matrix[::-1]))
    return pwm


def count(sequence, SC, ICPC, ML):
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
        counter[0] = -math.log2(counter[0])
    if counter[1] != 0.0:
        counter[1] = -math.log2(counter[1])
    if counter[2] != 0.0:
        counter[2] = -math.log2(counter[2])
    if counter[3] != 0.0:
        counter[3] = -math.log2(counter[3])

    sum_counter = sum(counter)
    if sum_counter == 0 or sum_counter == 0.0:
        for i in range(len(counter)):
            if counter[i] == -0.0:
                counter[i] = ICPC
    else:
        total = ICPC/sum_counter
        counter = [x * total for x in counter]

    return counter


def generate_binding_sites(pwm, SC, SL, ML):
    binding_sites = []
    for i in range(SC):
        binding_sites.append(generate_binding_site(pwm, SC, SL, ML))

    return binding_sites


def generate_binding_site(pwm, SC, SL, ML):
    bases = ['A', 'C', 'G', 'T']
    binding_site = ''

    for i in range(ML):
        random_y = random.randint(0, 3)
        random_x = random.randint(0, ML-1)

        num = math.floor((pwm[random_y][random_x] * random.randint(0, ML))) % 4
        if random.random() < .5:
            binding_site += bases[num]
    return binding_site


def plant(sequences, binding_sites, SL, ML):

    curr = 0
    planted_sites = []
    temp = []
    for row in sequences:
        rand_starting_index = random.randint(0, SL-ML)
        collapsed_row = ''.join(row)
        temp2 = copy.deepcopy(collapsed_row)

        #collapsed_row.replace(collapsed_row[rand_starting_index:(rand_starting_index+ML)], binding_sites[curr])
        collapsed_row = collapsed_row[:rand_starting_index] + binding_sites[curr] + collapsed_row[rand_starting_index+len(binding_sites[curr]):]
        planted_sites.append("{} {}".format(rand_starting_index, rand_starting_index + len(binding_sites[curr])))
        curr += 1
        temp.append(collapsed_row)



    new_sequences = []
    for row in temp:

        new_sequences.append(list(row))
    return new_sequences, planted_sites


def generate_data(ICPC, ML, SL, SC, dir):
    os.makedirs(dir, exist_ok=True)
    sequences = generate_sequence(SC, SL)
    pwm = generate_motif(SC, SL, ML, ICPC)
    binding_sites = generate_binding_sites(pwm, SC, SL, ML)
    planted, planted_sites = plant(sequences, binding_sites, SL, ML)

    # for j in range(len(sequences)):
    #     for i in range(len(sequences[0])):
    #         if sequences[j][i] != planted[j][i]:
    #             print("{} {} index and {} sequences and {} planted".format(j, i, sequences[j][i], planted[j][i]))

    with open(dir+"sequences.fa", "w+") as file:
        i = 1
        for row in planted:
            temp = ""
            for column in row:
                temp += column
            file.write(">Sequence {}\n".format(i))
            file.write(temp + "\n")
            i += 1

    with open(dir+"sites.txt", "w+") as file:
        i = 1
        for row in planted_sites:
            line = row.strip().split(" ")
            file.write("Start Index: {} End Index: {}\n".format(line[0], line[1]))

    zipped = zip(*pwm)

    with open(dir+"motif.txt", "w+") as file:
        i = 1
        file.write(">MOTIF{}\t{}\n".format(i, ML))
        for row in zipped:
            for column in row:
                file.write("{}\t".format(column))
            file.write("\n")
        file.write("<")

    with open(dir+"motiflength.txt", "w+") as file:
        file.write(str(ML))


def generate_sets(ICPC, ML, SL, SC, dir):
    os.makedirs(dir, exist_ok=True)
    for i in range(1, 11):
        generate_data(ICPC, ML, SL, SC, "{}/data{}/".format(dir, i))

# default values
# ICPC = 2
# ML = 8
# SL = 500
# SC = 10
# generate_data(2, 8, 500, 10, "generated_data/")
# generate_sets(2, 8, 500, 10, "generated_data/set1")

generate_sets(2, 8, 500, 10, "data_set/set1/")
generate_sets(1, 8, 500, 10, "data_set/set2/")
generate_sets(1.5, 8, 500, 10, "data_set/set3/")
generate_sets(2, 6, 500, 10, "data_set/set4/")
generate_sets(2, 7, 500, 10, "data_set/set5/")
generate_sets(2, 8, 500, 5, "data_set/set6/")
generate_sets(2, 8, 500, 20, "data_set/set7/")

import os
import math

def calc_overlap(predicted_sites_fname, sites_fname):
	predicted_sites_f = open(predicted_sites_fname, "r")
	sites_f = open(sites_fname, "r")

	predicted_sites_lines = predicted_sites_f.readlines()
	sites_lines = sites_f.readlines()
	overlap_cnt = 0

	for i in range(len(sites_lines)):
		site_line_split = sites_lines[i].split()
		predicted_site_line_split = predicted_sites_lines[i].split()

		site_start = int(site_line_split[2])
		site_end = int(site_line_split[5])

		predicted_start = int(predicted_site_line_split[2])
		predicted_end = int(predicted_site_line_split[5])

		if max(predicted_start, site_start) < min(site_end, predicted_end):
			overlap_cnt += 1

	return overlap_cnt

def get_background_probabilities(motif_matrix)
	q = {}
	a_cnt = 0
	c_cnt = 0
	g_cnt = 0
	t_cnt = 0
	for i in range(4):
		for j in range(len(motif_matrix)):
			if i == 0:
				a_cnt += motif_matrix[j][i]
			elif i == 1:
				c_cnt += motif_matrix[j][i]
			elif i == 2:
				g_cnt += motif_matrix[j][i]
			elif i == 3:
				t_cnt += motif_matrix[j][i]

	total_cnt = a_cnt + c_cnt + g_cnt + t_cnt
	q["A"] = float(a_cnt)/total 
	q["C"] = float(c_cnt)/total 
	q["T"] = float(t_cnt)/total 
	q["G"] = float(g_cnt)/total 

	return q

def create_pwm(M):
	pwm = [[0 for j in range(len(M[0]))] for i in range(len(M))]

	for i in range(len(M)):
		row_sum = 0
		for j in range(len(M[0])):
			row_sum += M[i][j]

		for j in range(len(M[0])):
			pwm[i][j] = float(M[i][j])/row_sum

	return pwm

def calc_entropy(motif_fname, predicted_motif_fname):
	motif_matrix = open(motif_fname, "r").readlines()
	motif_matrix = motif_matrix[1:len(motif_matrix)-1]
	for i in range(len(motif_matrix)):
		motif_matrix[i] = motif_matrix[i].strip()
		motif_matrix[i] = motif_matrix[i].split("\t")
		motif_matrix[i] = [int(x) for x in m[i]]

	q = get_background_probabilities(motif_matrix)

	predicted_matrix = open(predicted_motif_fname, "r").readlines()
	predicted_matrix = predicted_matrix[1:len(predicted_matrix)-1]
	for i in range(len(predicted_matrix)):
		predicted_matrix[i] = predicted_matrix[i].strip()
		predicted_matrix[i] = predicted_matrix[i].split("\t")
		predicted_matrix[i] = [int(x) for x in m[i]]

	W = create_pwm(predicted_matrix)

	entropy = 0
	for i in range(len(W)):
		for j in range(len(W[0])):
			base_prob = 0
			if j == 0:
				base_prob = q["A"]
			elif j == 1:
				base_prob = q["C"]
			elif j == 2:
				base_prob = q["G"]
			elif j == 3:
				base_prob = q["T"]

			entropy += (W[i][j]*math.log(W[i][j]/base_prob, 2))

	return entropy

overlap_cnt = 0
for i in range(1, 8):
	for j in range(1, 11):
		for folder, subs, files in os.walk("data_set\set"+str(i)+"\data"+str(j)):
			predicted_sites_fname = "data_set\set"+str(i)+"\data"+str(j)+"\predictedsites.txt"
			sites_fname = "data_set\set"+str(i)+"\data"+str(j)+"\sites.txt"
			overlap_cnt += calc_overlap(predicted_sites_fname, sites_fname)
			print "Number of overlapping sites for dataset"+str(i)+ ", data"+str(j)+": " + str(overlap_cnt)

			motif_fname = "data_set\set"+str(i)+"\data"+str(j)+"\motif.txt"
			predicted_motif_fname = "data_set\set"+str(i)+"\data"+str(j)+"\predictedmotif.txt"
			relative_entropy = calc_entropy(motif_fname, predicted_motif_fname)
			print "Entropy for dataset"+str(i)+ ", data"+str(j)+": " + str(relative_entropy)

			#Run time here
			






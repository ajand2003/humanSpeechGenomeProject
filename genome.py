

# Dictionary mapping Symbol to Name
amino_acids = {"A": "Alanine", "I": "Isoleucine",
	"L": "Leucine", "V": "Valine", "F": "Phenylalanine", 
	"W": "Tryptophan", "Y": "Tyrosine", "N": "Asparagine",
	"C": "Cystenine", "Q": "Glutamine", "T": "Theronine", 
	"D": "Aspartic Acid", "E": "Glutamic Acid", "R": "Arginine",
	"H": "Histidine", "K": "Lysine", "G": "Glycine", "P": "Proline"}

def parse_genome(file) -> []:
	file = open(file, "r")
	genome_list = []
	for line in file:
		sequence = line.replace(" ", "").replace("\n", "")
		for amino_acid in sequence:
			genome_list.append(amino_acid)
	print(len(genome_list))
	return genome_list


def lcs(seq1,seq2) -> []:
	len1 = len(seq1)
	len2 = len(seq2)

	matrix = [[0] * (len2 + 1) for _ in range(len1 + 1)]

	for i in range(1, len1 + 1):
		for j in range(1, len2 + 1):
			if seq1[i-1] == seq2[j-1]:
				matrix[i][j] = matrix[i-1][j-1] + 1
			else:
				matrix[i][j] = max(matrix[i][j-1], matrix[i-1][j])
	return matrix 
		
def get_lcs(C, seq1, seq2):
    result = ""
    i = len(seq1)
    j = len(seq2)

    while i > 0 and j > 0:
    	if C[i][j] == C[i-1][j]:
    		i = i - 1
    	elif C[i][j] == C[i][j-1]:
    		j = j - 1
    	else:
    		assert seq1[i-1] == seq2[j-1]
    		result = seq1[i-1] + result
    		i = i - 1
    		j = j - 1
    return result

def find_differences(C, seq1, seq2):
	result = ""
	i = len(seq1)
	j = len(seq2)

	while i > 0 and j > 0:
		if seq1[i-1] == seq2[j-1]:
			result = seq1[i-1] + result 
			i = i - 1
			j = j - 1
		else:
			if C[i][j-1] >= C[i-1][j]:
				result = " +" + seq2[j-1] + " " + result
				j = j - 1
			elif C[i][j-1] < C[i-1][j]:
				result = " -" + seq1[i-1] + " " + result
				i = i - 1
	return result

def compare_genomes(seq1,seq2) -> []: 
	index1 = 0
	index2 = 0 
	lst = lcs(seq1, seq2)
	return lst

def number_of_differences(seq):
	positive_count = 0
	negative_count = 0
	differences = 0
	print(seq)
	for i in range(len(seq)):
		if seq[i] == "+":
			positive_count = positive_count + 1
		elif seq[i] == "-":
			negative_count = negative_count + 1
	if positive_count > negative_count:
		differences = positive_count
	elif negative_count > positive_count:
		differences = negative_count
	else:
		differences = positive_count
	return differences

def print_table(C):
	dictionary = {0: "gibbon", 1: "gorilla", 2: "human", 3: "housemouse",
				  4: "bonobo", 5: "commonchimp", 6: "orangutan",
				  7: "zebrafish", 8: "hummingbird", 9: "americanalligator",
				  10: "westernclawedfrog"}
	for i in range(len(C)):
		for j in range(len(C[i])):
			print(str(C[i][j]) + " " + 3 * " ", end="")
		print("\n")

def generate_tree(C):
	dictionary = {0: "gibbon", 1: "gorilla", 2: "human", 3: "housemouse",
				  4: "bonobo", 5: "commonchimp", 6: "orangutan",
				  7: "zebrafish", 8: "hummingbird", 9: "americanalligator",
				  10: "westernclawedfrog"}
	C[0][0] = 1000
	min_pos = (0,0)
	while input() != "q":
		print_table(C)
		for i in range(len(C)):
			for j in range(len(C[i])):
				if (C[i][j] < C[min_pos[0]][min_pos[1]]):
					min_pos = (i,j)
		C = calculate_new_table(C, min_pos[0], min_pos[1])
	print(min_pos)

def average(numbers):
	return float(sum(numbers)) / max(len(numbers), 1)

def calculate_new_table(C, group1, group2):
	new_table = []
	lst = []
	for j in range(len(C[group1])):
		lst.append(average([C[group1][j], C[group2][j]]))
	new_table.append(lst)
	for i in range(len(C)):
		if not group1 == i and not group2 == i:
			lst = []
			for j in range(len(C[i])):
				lst.append(C[i][j])
			new_table.append(lst)
	return new_table

def get_differences(lists):
	dictionary = {0: "gibbon", 1: "gorilla", 2: "human", 3: "housemouse",
				  4: "bonobo", 5: "commonchimp", 6: "orangutan",
				  7: "zebrafish", 8: "hummingbird", 9: "americanalligator",
				  10: "westernclawedfrog", 11: "dog"}
	C = []
	for i in range(len(lists)):
		row_list = []
		for j in range(len(lists)):
			lst = compare_genomes(lists[i], lists[j])
			common_sequence = find_differences(lst, lists[i], lists[j])
			differences = number_of_differences(common_sequence)
			print(dictionary[i] + " ---> " + dictionary[j] + " " + str(differences))
			if j > i:
				row_list.append(differences)
			else:
				row_list.append(0)
		C.append(row_list)
	generate_tree(C)
	print(C)


#Main Entry point
if __name__ == "__main__":
	gibbon_seq = parse_genome("realgibbonfoxp2.txt")
	gorilla_seq = parse_genome("realgorillafoxp2.txt")
	human_seq = parse_genome("realhomosapienfoxp2.txt")
	housemouse_seq = parse_genome("realhousemousefoxp2.txt") 
	bonobo_seq = parse_genome("realpanpaniscusfoxp2.txt")
	commonchimp_seq = parse_genome("realpantrogladytefoxp2.txt")
	orangutan_seq = parse_genome("realpongofoxp2.txt")
	zebrafish_seq = parse_genome("realzebrafishfoxp2.txt")
	hummingbird_seq = parse_genome("realhummingfoxp2.txt")
	americanalligator_seq = parse_genome("realamericanalligatorfoxp2.txt")
	westernclawedfrog_seq = parse_genome("realwesternclawedfrogfoxp2.txt")
	dog_seq = parse_genome("realdogfoxp2.txt")
	'''lst = []
	lst.append(gibbon_seq)
	lst.append(gorilla_seq)
	lst.append(human_seq)
	lst.append(housemouse_seq)
	lst.append(bonobo_seq)
	lst.append(commonchimp_seq)
	lst.append(orangutan_seq)
	lst.append(zebrafish_seq)
	lst.append(hummingbird_seq)
	lst.append(americanalligator_seq)
	lst.append(westernclawedfrog_seq)
	lst.append(dog_seq)
	#get_differences(lst)
	C = compare_genomes(zebrafish_seq, hummingbird_seq)
	#print(get_lcs(C, gibbon_seq, hummingbird_seq))
	#print("\n")
	#print(find_differences(C,gibbon_seq, hummingbird_seq))'''
 


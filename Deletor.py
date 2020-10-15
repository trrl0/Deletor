import os
import itertools


def generate_sequences(start_sequence):
# Generates a list of all possible sequences with different amino acids deleted
    master_list = [] # Hold list of lists with individual elements of residue letters
    residues_sequence = start_sequence[:-1] # Leaves off the last amino acid which was the resin
    resin = start_sequence[-1]

    aloop = 1 # Starting with all combinations of 1 number
    bloop = 0 # Position in master_list for adding to current list by extracting values from the tuples in itertools

    while aloop <= len(residues_sequence): # Extra list manipulation is due to itertools returning lists of tuples
        for item in list(itertools.combinations(residues_sequence,aloop)):
            master_list.append([])
            for tuple_pos in item:
                master_list[bloop].append(tuple_pos)
            bloop += 1
        aloop += 1

    sequence_str = ""
    sequence_list = []
    for combo_list in master_list:
        for residue in combo_list:
            sequence_str += residue
        sequence_str += resin
        sequence_list.append(sequence_str)
        sequence_str = ""

    return sequence_list



def generate_masses(products):

    amino_masslist = {
        'A': 71.03711,
        'R': 156.10111,
        'N': 114.04293,
        'D': 115.02694,
        'C': 103.00919,
        'E': 129.04259,
        'Q': 128.05858,
        'G': 57.02146,
        'H': 137.05891,
        'I': 113.08406,
        'L': 113.08406,
        'K': 128.09496,
        'M': 131.04049,
        'F': 147.06841,
        'P': 97.05276,
        'S': 87.03203,
        'T': 101.04768,
        'W': 186.07931,
        'Y': 163.06333,
        'V': 99.06841,

        # Special residues

        'U': (129.04259 + 248.05),  # Edans
        'B': (128.09496 + 251.11),  # Dabcyl
        'O': (128.05858 - 17.03),   # pyroGlu
        'X': (133.06),              # N-MeAbz on N-term
        'Z': (165.01 + 87.03203),   # Dinitrophenol
        '4': 229.93                 # 4-iodobenzoic acid

    }

    mass = 0
    mass2 = 0  # 2+ charge state
    mass3 = 0
    mass4 = 0
    peptide_masses = []
    for peptide in products:
        for amino in peptide:
            mass += amino_masslist[amino]
        mass += 19.01784  # N-term, C-term, proton
        mass2 = (mass + 1.00783) / 2  # Second proton
        mass3 = (mass + 2 * 1.00783) / 3
        mass4 = (mass + 3 * 1.00783) / 4
        peptide_masses.append([peptide, mass, mass2, mass3, mass4])
        mass = 0
        mass2 = 0

    return peptide_masses

def save_masslist(peptide_masses, filename):
    filename += ".csv"
    outStr = ""

    try:
        print ("Saving " + filename + "...", end="", flush=True)

        outStr += "Sequence,[M+H]+,[M+2H]2+,[M+3H]3+,[M+4H]4+\n"  # Headers for excel

        for item in peptide_masses:
            outStr += str(item[0]) + "," + str(item[1]) + "," + str(item[2]) + "," + str(item[3]) + "," + str(
                item[4]) + "\n"

        out_file = open(os.getcwd() + "\\" + filename, "w")
        out_file.write(outStr)
        out_file.close()

        print ("done")
    except:
        print("File saving error! File may be open or otherwise locked by the OS.")


parent_sequence = input("Enter the peptide sequence: ")
save_name = parent_sequence[0:3] + "_deletelist"


save_masslist(
     generate_masses(
        generate_sequences(parent_sequence)
     ),
       save_name
)

input ("Press ENTER to exit...")


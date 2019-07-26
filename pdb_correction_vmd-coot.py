#!/usr/bin/env python3

import argparse
import ipdb

def number_to_hybrid36_number( number, width):
		
	digits_upper = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"

	def encode_pure(digits, value):
		"encodes value using the given digits"
		assert value >= 0
		if (value == 0): return digits[0]
		n = len(digits)
		result = []
		while (value != 0):
			rest = value // n
			result.append(digits[value - rest * n])
			value = rest
		result.reverse()
		return "".join(result)

	#digits_upper_values = dict([pair for pair in zip(digits_upper, range(36))])
	digits_lower = digits_upper.lower()
	#digits_lower_values = dict([pair for pair in zip(digits_lower, range(36))])
	if (number >= 1-10**(width-1)):
		if (number < 10**width):
			return '{:{width}d}'.format(number, width=width)
		number -= 10**width
		if (number < 26*36**(width-1)):
			number += 10*36**(width-1)
			return encode_pure(digits_upper, number)
		number -= 26*36**(width-1)
		if (number < 26*36**(width-1)):
			number += 10*36**(width-1)
			return encode_pure(digits_lower, number)
	raise ValueError("value out of range.")

REPLACEMENT_DICT = {' O1P':' OP1', ' O2P':' OP2', ' C5M':' C7 '}
REPLACEMENT_DICT_BASES = {'CYT':' DC', 'GUA':' DG', 'THY':' DT', 'ADE':' DA', 'DA5':' DA', 'DA3':' DA', 'DT5':' DT', 'DT3':' DT', 'DG5':' DG', 'DG3':' DG', 'DC5':' DC', 'DC3':' DC'}
OCC = "9.99"
BFAC= "1.00"

#TODO: -high cleanup
class PDB_Corr(object):
	
	def __init__(self):
		self.current = {"atom_number": 1, "molecule_number": 1, "last_molecule_number": 1, "chain_id": "A", "chain_id_repeats": 1, "chain": None}

	def reshuffle_pdb(self, pdb_file):
		unshuff_file = []
		for line in pdb_file:
			unshuff_file.append(line)
		####
		shuff_file = unshuff_file
		####
		newFile = ''.join(shuff_file)
		return newFile

	def correct_pdb(self, pdb_file):
		reset_numbers = True
		header = "AUTHORS:     Martin,Casanal,Feigl\n"
		#TODO: -low MODEL ?
		body = [header]
		for line in pdb_file:
			lineType = line[0:6]
			if lineType in ["TITLE ", "CRYST1"]:
				body.append(line)
			if lineType == "ATOM  ":
				no_atom_check = line[13:15]
				if not(no_atom_check == "  "):
					#TODO: -mid restore options
					line = self.correct_atom_number(line)
					line = self.correct_nomenclature(line)
					line, is_ter = self.correct_molecule_chain_and_number(line, reset_numbers)
					line = self.correct_occupancy(line)
					line = self.correct_atomtype(line)
					if is_ter:
						body.append("TER\n")
					body.append(line)
		body.append("END\n")
		newFile = ''.join(body)
		return newFile



	def correct_atom_number(self, line):
		atom_number_string = number_to_hybrid36_number(self.current["atom_number"], 5)
		newline = line[0:5] + " " + atom_number_string + line[11:]
		self.current["atom_number"] += 1
		return newline

	def correct_nomenclature(self, line):
		atom = line[12:16]
		if atom in REPLACEMENT_DICT:
			atom = REPLACEMENT_DICT.get(atom, '    ')
		base = line[17:20]
		if base in REPLACEMENT_DICT_BASES:
			base = REPLACEMENT_DICT_BASES.get(base, '   ')
		newline = line[0:12] + atom + ' ' + base + line[20:]
		return newline
			
	def correct_occupancy(self, line):
		newline = line[:54] + "  " + BFAC + "  " + OCC + line[66:]
		return newline

	def correct_molecule_chain_and_number (self, line, reset_numbers):
		def increase_chain_id(current_chain_id):
			possible_chain_ids = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
			pos = possible_chain_ids.find(current_chain_id)
			chlen = len(possible_chain_ids)
			if (pos +1 < chlen):
				return possible_chain_ids[pos+1]
			else:
				# Don't allow for A again (Scaffold) that's why 1 and not 0
				return possible_chain_ids[1]

		chain_id = line[21:22]
		chain = line[72:76]
		is_ter = False
		
		molecule_number_str = line[22:26]
		molecule_number = int(molecule_number_str.replace(' ',''))
		
		chain_id_only = True if chain == "    " else False

		if chain_id_only: #TODO
			if chain_id != self.current["chain_id"]:
				new_chain_id = increase_chain_id(chain_id)
				new_molecule_number = 1
							
			elif molecule_number != self.current["molecule_number"]:
				new_molecule_number += 1
			
			else:
				new_chain_id = self.current["chain_id"]
				new_molecule_number += 1

		else:
			if self.current["chain"] is None:
				self.current["chain"] = chain
				new_chain_id = self.current["chain_id"]
				new_molecule_number = 1
			elif chain == self.current["chain"]:
				new_chain_id = self.current["chain_id"]  
				if molecule_number == self.current["molecule_number"]:
					new_molecule_number = self.current["last_molecule_number"]
				else:
					new_molecule_number = self.current["last_molecule_number"] +1
			else:
				new_chain_id = increase_chain_id(self.current["chain_id"])
				is_ter = True
				new_molecule_number = 1
				if self.current["chain_id"] == "Z":
					self.current["chain_id_repeats"] += 1

		new_chain_str = str(new_chain_id) + str(self.current["chain_id_repeats"]).rjust(3,"0")
		new_molecule_number_str = number_to_hybrid36_number(1000*(self.current["chain_id_repeats"]-1) + new_molecule_number, 4)
		newline = line[0:21] + new_chain_id + new_molecule_number_str + line[26:67]+ "     " + new_chain_str + line[76:]

		self.current["chain_id"] = new_chain_id
		self.current["chain"] = chain
		self.current["molecule_number"] = molecule_number
		self.current["last_molecule_number"] = new_molecule_number
		
		return newline, is_ter

	def correct_atomtype(self, line):
		atom = line[12:14]
		if atom[1] in "0123456789":
			atom = atom[0]

		atom = atom.strip().rjust(2," ")
		newline = line[:76] + atom + line[78:]
		#print(atom)
		return newline


#TODO: -low restore remove H funktionality

def main():
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('--i',nargs='?', help='input filename',type=str)
	parser.add_argument('--o',nargs='?', help='output filename',type=str)
	parser.add_argument('--h', help='keep header',action="store_true")
	parser.add_argument('--f', help='keep footer',action="store_true")
	#TODO: -mid restore flexibility

	args = parser.parse_args()

	print(args.i)
	print ("start")
	pdb_Corr = PDB_Corr()
	with open(args.i,'r') as file_init:
		reshuffFile = pdb_Corr.reshuffle_pdb(file_init)

	with open(args.o, "w") as file_reshuff:
		file_reshuff.write(reshuffFile)

	with open(args.o, "r") as file_reshuff:
		newFile = pdb_Corr.correct_pdb(file_reshuff)

	with open(args.o, "w") as file_corr:
		file_corr.write(newFile)
		
	print ("done")
	return


if __name__ == "__main__":
    main()

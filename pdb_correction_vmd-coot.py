#!/usr/bin/env python
# coding: utf-8

import argparse

version = "0.1.0"
header = "HEADER\n AUTHORS		: Martin,Casanal\n"


def correct_pdb (pdb_file, keep_header, keep_footer, nomenclature, remove_h, molecule_chain, atom_number, occupancy):
	started = False
	newFileArray = []
	footer = []
	if not (keep_header):
		newFileArray.append(header)
	for line in pdb_file:
		lineType = line[0:4]
		if (lineType == "ATOM"):
			started = True
			footer = []
			no_atom_check = line[13:15]
			if not(no_atom_check == "  "):
				# correct line
				append = True

				if (nomenclature):
					line = correct_nomenclature(line)
				if (remove_h):
					if not correct_remove_H_if_False(line):
						append = False
				if (molecule_chain):
					line = correct_molecule_chain_and_number(line)
				if (atom_number):
					line = correct_atom_number(line)
				if (occupancy):
					correct_occupancy

				if (append):
					newFileArray.append(line)
			
		elif not started:
			if (keep_header):
				newFileArray.append(line)
		else:
			footer.append(line)
	if not (keep_footer):
		newFileArray.append('END')
	else:
		newFileArray.extend(footer)

	newFile = ''.join(newFileArray)
	return newFile

atom_number = 1
molecule_number = 0
chain = "A"
last_molecule_number = 0
last_chain_id = ""
chain_id_repeats = 0
def reset_global_variables ():
	global atom_number
	global molecule_number
	global chain
	global last_molecule_number
	global last_chain_id
	global chain_id_repeats
	atom_number = 1
	molecule_number = 0
	chain = "A"
	last_molecule_number = 0
	last_chain_id = ""
	chain_id_repeats = 0


def correct_atom_number (line):
	global atom_number
	atom_number_string = number_to_hybrid36_number(atom_number, 5)
	newline = line[0:5] + " " + atom_number_string + line[11:]
	atom_number += 1
	return newline


digits_upper = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
digits_upper_values = dict([pair for pair in zip(digits_upper, range(36))])
digits_lower = digits_upper.lower()
digits_lower_values = dict([pair for pair in zip(digits_lower, range(36))])
def number_to_hybrid36_number(number, width):
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


def correct_molecule_chain_and_number (line):
	global molecule_number
	global chain
	global last_molecule_number
	global last_chain_id
	global chain_id_repeats

	current_chain_id = line[21:22]
	current_molecule_number_str = line[22:26]
	current_molecule_number = int (current_molecule_number_str.replace(' ',''))
	if (current_molecule_number < last_molecule_number):
		chain = increase_chain_id(chain)
		molecule_number = 1
	elif not (current_molecule_number == last_molecule_number):
		molecule_number += 1
  
	if (last_chain_id == "Z") and not (last_chain_id == chain):
		chain_id_repeats += 1

	chain_repeat_string = " "
	if (chain_id_repeats > 0):
		chain_repeat_string = str(chain_id_repeats)

	molecule_string = number_to_hybrid36_number(1000*chain_id_repeats + molecule_number, 4)
	newline = line[0:21] + chain + molecule_string + line[26:72] + chain + chain_repeat_string + "  " + line[76:]

	last_molecule_number = current_molecule_number
	last_chain_id = chain

	return newline
			

possible_chain_ids = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
def increase_chain_id (current_chain_id):
	pos = possible_chain_ids.find(current_chain_id)
	chlen = len(possible_chain_ids)
	if (pos +1 < chlen):
		return possible_chain_ids[pos+1]
	else:
		# Don't allow for A again (Scaffold) thats why 1 and not 0
		return possible_chain_ids[1]

#def number_to_ndigit_string (number, width):
#	return '{:{width}d}'.format(number, width=width)


replacement_dict = {' OP1':' O1P', ' OP2':' O2P'}
replacement_dict_bases = {'CYT':' DC', 'GUA':' DG', 'THY':' DT', 'ADE':' DA', 'DA5':' DA', 'DA3':' DA', 'DT5':' DT', 'DT3':' DT', 'DG5':' DG', 'DG3':' DG', 'DC5':' DC', 'DC3':' DC'}
def correct_nomenclature (line):
	atom = line[12:16]
	if atom in replacement_dict:
		atom = replacement_dict.get(atom, '    ')
	base = line[17:20]
	if base in replacement_dict_bases:
		base = replacement_dict_bases.get(base, '   ')
	newline = line[0:12] + atom + ' ' + base + line[20:]
	return newline
			

def correct_remove_H_if_False (line):
	atom = line[12:16]
	if not(atom == "    ") and "H" not in atom:
		return True
	return False


def correct_occupancy (line):
	newline = line[:54] + '  1.00  0.00' + line[67:]
	return newline
				

parser = argparse.ArgumentParser(description='')

parser.add_argument('--i',nargs='?', help='input filename',type=str)
parser.add_argument('--o',nargs='?', help='output filename',type=str)
parser.add_argument('--h', help='keep header',action="store_true")
parser.add_argument('--f', help='keep footer',action="store_true")


args = parser.parse_args()

print(args.i)

print ("start")

file=open(args.i,'r')
newFile = correct_pdb(file,args.h,args.f,True,True,True,True,True)
text_file = open(args.o, "w")
text_file.write(newFile)
text_file.close()

print ("done")





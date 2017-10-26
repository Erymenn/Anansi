#!/usr/bin/python -i

import os
import sys
import networkx as nx
import numpy as np


def data_to_dump(datafile, dumpfile, zones = None, dumpbonds = None):
	"""(str, str, [int], str) -> None
	Write a dump Lammps file from a data Lammps file with columns id mol type x y z
	If zones (list of units id), rewrite the mol numbers to match zones
	If dumpbonds, use the dumpbonds file to mark the broken bond by changing the type of their units
	"""
	graph = data_to_graph(datafile)
	if dumpbonds:
		brokenbonds = dumpbonds_to_brokenbondssets(dumpbonds)[-1]
		for a, b in brokenbonds:
			graph.node[a]['type'] = INDEX['broken']
			graph.node[b]['type'] = INDEX['broken']
	if zones:
		#creation of dictionnary id-mol
		dico = {z:i for i,z in enumerate(zones, 2)}
		#mol replacement
		for node in graph.nodes_iter():
			if node in dico: graph.node[node]['mol'] = dico[node]
			else: graph.node[node]['mol'] = 1
	graph_to_dump(graph, dumpfile)


def data_to_graph(filename):
	"""(str) -> nx.Graph
	Return the graph based on the filname data Lammps file
	This graph include:
		attributs: heading, masses, pair and bonds coeffs, atom and bond types, atoms and bonds, extra bond per atom, xlo, xhi...
		nodes: units id with attributs mol, type, position, nodeimage, mask and velocity
		edges: bond with attribut weight which stands for the bond type
	"""
	#Creation of graph
	graphresult = nx.Graph()
	with open(filename,'r') as f:
		#Save of heading lines as graph attribut
		graphresult.graph['heading'] = f.readline().rstrip()
		graphresult.graph['masses'] = []
		graphresult.graph['pair coeffs'] = []
		graphresult.graph['bond coeffs'] = []
		f.readline()
		line = 1
		while line:
			line = f.readline().rstrip()
			if line: 
				line = line.split()
				graphresult.graph[" ".join(line[1:])] = int(line[0])
		for i in xrange(3):
			line = f.readline().rstrip().split()
			graphresult.graph[line[2]] = float(line[0])
			graphresult.graph[line[3]] = float(line[1])
		for i in xrange(3): f.readline()
		for i in xrange(graphresult.graph['atom types']):
			graphresult.graph['masses'].append(f.readline().rstrip())
		for i in xrange(3): f.readline()
		for i in xrange(graphresult.graph['atom types']):
			graphresult.graph['pair coeffs'].append(f.readline().rstrip())
		for i in xrange(3): f.readline()
		for i in xrange(graphresult.graph['bond types']):
			graphresult.graph['bond coeffs'].append(f.readline().rstrip())
		for i in xrange(3): f.readline()
		#Save the atoms tag as graph nodes
		#and atoms properties as graph attributs
		for i in xrange(graphresult.graph['atoms']):
			line = f.readline().rstrip().split()
			if i == graphresult.graph['atoms']: print line
			tag = int(line[0])
			graphresult.add_node(tag, mol = int(line[1]), type = int(line[2]))
			graphresult.node[tag]['position'] = np.array([float(line[3]), float(line[4]), float(line[5])])
			graphresult.node[tag]['nodeimage'] = " ".join(line[6:])
			graphresult.node[tag]['mask'] = False
		for i in xrange(3): f.readline()
		for i in xrange(graphresult.graph['atoms']):
			line = f.readline().rstrip().split()
			graphresult.node[int(line[0])]['velocity'] = np.array([float(line[1]), float(line[2]), float(line[3])])
		for i in xrange(3): f.readline()
		#Save the bonds as graph edges
		#and bond type as edge weight
		for i in xrange(graphresult.graph['bonds']):
			line = f.readline().rstrip().split()
			graphresult.add_edge(int(line[2]), int(line[3]), type = int(line[1]))
	return graphresult


def dump_to_sequencegraphswithoutbonds(filename):
	"""(str) -> [nx.Graph]
	Return a list of graph built each dump timestep
	dump file columns: id type mol x y z
	each graph include:
		attributs: timestep and xlo, xhi, ylo...
		nodes: units id with attributes type, mol et position
		no edges
	"""
	sequencegraphswithoutbonds = []
	with open(filename,'r') as f:
		while f.readline():
			#as long as something can be read
			g = nx.Graph()
			g.graph["timestep"] = f.readline().rstrip()
			f.readline()
			nbatt = int(f.readline().rstrip())
			f.readline()
			for n in ('x', 'y', 'z'):
				g.graph[n+'lo'], g.graph[n+'hi'] = [float(i) for i in f.readline().rstrip().split()]
			cols = f.readline().split()
			for i in xrange(nbatt):
				line = f.readline().split()
				tag = int(line[0])
				g.add_node(tag, type = int(line[1]), mol = int(line[2]))
				g.node[tag]['position'] = np.array([float(line[3]), float(line[4]), float(line[5])])
			sequencegraphswithoutbonds.append(g)
	return sequencegraphswithoutbonds


def dumpbonds_to_brokenbondssets(dump):
	"""(str) -> [set()]
	Return a list of sets containing the bonds broen scince the last timestep
	3 first columns of dump file: id1 id2 type
	"""
	sets = list()
	precnbbd = -1
	bds = set()
	with open(dump,'r') as f:
		while f.readline():
			f.readline()
			f.readline()
			nbbd = int(f.readline().rstrip())
			for i in xrange(5): f.readline()
			if nbbd == precnbbd:
				sets.append(set())
				for i in xrange(nbbd): f.readline()
			else:
				newbds = set()
				if precnbbd == -1:
					sets.append(set())
					for i in xrange(nbbd): 
						l = f.readline().split()
						bds.add((int(l[0]), int(l[1])))
				else:
					for i in xrange(nbbd): 
						l = f.readline().split()
						newbds.add((int(l[0]), int(l[1])))
					sets.append(bds.difference(newbds))
					bds = newbds
					assert len(sets[-1]) == precnbbd-nbbd
				precnbbd = nbbd
	return sets


def dumpbonds_to_bondssets(dump):
	"""(str) -> [set()]
	Return a list including the bond for each dump timestep
	3 first columns of dump file: id1 id2 type
	"""
	sets = list()
	with open(dump,'r') as f:
		while f.readline():
			s = set()
			f.readline()
			f.readline()
			nbbd = int(f.readline().rstrip())
			for i in xrange(5): f.readline()
			for i in xrange(nbbd):
				l = f.readline().split()
				s.add((int(l[0]), int(l[1])))
			sets.append(s)
	return sets
	
	
def dumpbonds_to_bondssets2(dump):
	"""(str) -> [set()]
	Return a list including the bond for every dump timestep
	3 first columns of dump file: id1 id2 type
	"""
	sets = list()
	s = set()
	with open(dump,'r') as f:
		for l in f:
			if l.startswith("ITEM: ATOMS"): 
				lirebd = True
			elif lirebd and l.startswith('ITEM'): 
				lirebd = False
				sets.append(s)
				s = set()
			elif lirebd:
				l = l.split()
				s.add((int(l[0]), int(l[1])))
		sets.append(s)
	return sets


def dumpbonds_to_weighted_list(dump):
	"""(str, bool) -> [list()]
	Return a list of built bonds for every dump timestep
	columns of dump file: id1 id2 type distance energy force
	"""
	wl = list()
	truc = 0
	with open(dump,'r') as f:
		while f.readline():
			s = list()
			f.readline()
			f.readline()
			nbbd = int(f.readline().rstrip())
			for i in xrange(5): f.readline()
			for i in xrange(nbbd):
				l = f.readline().split()
				s.append((int(l[0]), int(l[1]), int(l[2]), float(l[3]), float(l[4]), float(l[5])))
			wl.append(s)
	return wl
	

def dumps_to_sequencegraphwithbonds(dump, dumpbonds):
	"""(str, str) -> [nx.Graph]
	Return a list of graph built for every dump timestep
	columns of dump file: id type molecule x y z
	columns of dumpbonds file: id1 id2 type distance energy force 
	Each graph include:
		attributs: timestep and xlo, xhi, ylo...
		nodes: units tag with attributs type, mol and position
		edges: bonds with attributs distance, energy and force
	"""
	sequencegraphwithbonds = dump_to_sequencegraphswithoutbonds(dump)
	edges = dumpbonds_to_weighted_list(dumpbonds)
	for n in xrange(len(sequencegraphwithbonds)): sequencegraphwithbonds[n].add_weighted_edges_from(edges[n])
	return sequencegraphwithbonds


def graph_to_data(graph, name):
	"""(nx.Graph, str) -> None
	Write a data file in current directory
	The file is build from graph and its name is name
	"""
	if not name in os.listdir(os.getcwd()): open(name, 'w').close()
	if 'extra bond per atom' not in graph.graph: graph.graph['extra bond per atom'] = 4
	with open(name,'w') as f:
		#print "...ecriture de l'en-tete"
		f.write(graph.graph['heading']+"\n\n")
		for att in graph.graph:
			if isinstance(graph.graph[att], int):
				f.write(str(graph.graph[att])+" "+att+"\n")
		f.write("\n")
		for i in ('x', 'y', 'z'):
			f.write("%.16e %.16e %slo %shi\n" %(graph.graph[i+'lo'], graph.graph[i+'hi'], i, i))
		f.write("\nMasses\n\n")
		f.write("\n".join(graph.graph["masses"]))
		f.write("\n\nPair Coeffs\n\n")
		f.write("\n".join(graph.graph["pair coeffs"]))
		f.write("\n\nBond Coeffs\n\n")
		f.write("\n".join(graph.graph["bond coeffs"]))
		f.write("\n\nAtoms\n\n")
		for n in graph.nodes_iter():
			f.write("%i %i %i %.16e %.16e %.16e %s\n" %(n, graph.node[n]["mol"], graph.node[n]["type"], graph.node[n]["position"][0], graph.node[n]["position"][1], graph.node[n]["position"][2], graph.node[n]["nodeimage"]))
		f.write("\nVelocities\n\n")
		for n in graph.nodes_iter():
			f.write("%i %.16e %.16e %.16e\n" %(n, graph.node[n]["velocity"][0], graph.node[n]["velocity"][1], graph.node[n]["velocity"][2]))
		f.write("\nBonds\n\n")
		i = 1
		for n1, n2 in graph.edges_iter():
			f.write("%i %i %i %i\n" %(i, graph.edge[n1][n2]['type'], n1, n2))
			i += 1


def graph_to_dump(graph, name):
	"""(nx.Graph, str) -> Graph
	Write a dump file in current directory
	The file is build from graph and its name is name
	colums of dump file: id type mol x y z   
	"""
	if not name in os.listdir(os.getcwd()): open(name, 'w').close()
	with open(name,'w') as f:
		f.write("ITEM: TIMESTEP\n")
		if "timestep" in graph.graph: f.write(graph.graph['timestep']+"\n")
		else: f.write("0\n")
		f.write("ITEM: NUMBER OF ATOMS\n"+str(len(graph.nodes())))
		f.write("\nITEM: BOX BOUNDS pp pp pp\n")
		for i in ('x', 'y', 'z'):
			f.write("%.5f %.5f\n" %(graph.graph[i+'lo'], graph.graph[i+'hi']))
		f.write("ITEM: ATOMS id type mol x y z\n")
		for n in graph.nodes_iter():
			f.write("%i %i %i %.5f %.5f %.5f\n" %(n, graph.node[n]["type"], graph.node[n]["mol"], graph.node[n]["position"][0], graph.node[n]["position"][1], graph.node[n]["position"][2]))
	
	
def renumber_id_data(filename, outputname = None):
	"""(str, str) -> None
	Write a data file from the file filename with the units id from 1 to (units number)
	outputname is the output file name, if no precise, the input file is rewritten
	"""
	g = data_to_graph(filename)
	dico = {k: n for n, k in enumerate(g.nodes_iter(),1)}
	g = nx.relabel_nodes(g, dico)
	if outputname: graph_to_data(g, outputname)
	else: graph_to_data(g, filename)


def renumber_id_dump(filename, outputname = None):
	"""(str, str) -> None
	Write a dump file from the file filename with the units id from 1 to (units number)
	outputname is the output file name, if no precise, the input file is rewritten
	"""
	seq = dump_to_sequencegraphswithoutbonds(filename)
	for n in xrange(len(seq)):
		dico = {k: i for i, k in enumerate(seq[n].nodes_iter(),1)}
		seq[n] = nx.relabel_nodes(seq[n], dico)
	if outputname: sequencegraphswithoutbonds_to_dump(seq, outputname)
	else: sequencegraphs_to_dump(sequencegraphswithoutbond, filename)
 
 
def sequencegraphs_to_dump(sequencegraphs, name):
	"""([nx.Graph], str) -> None
	Return a dump file named name for a list of graphs
	columns of dump: id type mol x y z
	Each graph include at least:
		attributs: timestep et xlo, xhi, ylo...
		nodes: units tag with attributs type, mol and position
	"""
	if not name in os.listdir(os.getcwd()): open(name, 'w').close()
	with open(name,'w') as f:
		for g in sequencegraphs:
			f.write("ITEM: TIMESTEP\n")
			f.write(g.graph['timestep']+"\n")
			f.write("ITEM: NUMBER OF ATOMS\n"+str(g.number_of_nodes()))
			f.write("\nITEM: BOX BOUNDS pp pp pp\n")
			for i in ('x', 'y', 'z'):
				f.write("%.5f %.5f\n" %(g.graph[i+'lo'], g.graph[i+'hi']))
			f.write("ITEM: ATOMS id type mol x y z\n")
			for n in g.nodes_iter():
				f.write("%i %i %i %.5f %.5f %.5f\n" %(n, g.node[n]["type"], g.node[n]["mol"], g.node[n]["position"][0], g.node[n]["position"][1], g.node[n]["position"][2]))


def slurm_to_follow_crosslinking(filename, outputname, initmonomers = 128000):
	"""(str, str, int, int) -> None
	Write a file named outputname which gave the progress of RLP polymerisation
	initmonomers is the initial number of monomers before RLP
	"""
	with open(filename,'r') as f:
		progressRLP = [initmonomers]
		next = False
		line = True
		while line:
			line = f.readline()
			if line.startswith('Molecule'): 
				next = True
			elif next:
				next = False
				line = line[:15].split()
				res.append(int(line[1][:-1]))
	t = [300*i for i in xrange(len(res))]
	print res, t
	finalnbbd = res[0]
	if not outputname in os.listdir(os.getcwd()): open(outputname, 'w').close()
	with open(outputname,'w') as f:
		f.write("#time(s) nbmonomers percentcrosslink\n")
		for i in xrange(len(t)):
			f.write('%d %d %f\n' %(t[i], res[i], float(res[0]-res[i])/finalnbbd))

#debug_graph = data_to_graph("../lib/graph/debug_data")
	
if __name__ == '__main__':
	#test
	debug_graph = data_to_graph("debug_data")
	graph_to_dump(debug_graph, "debug_dump_output")
	graph_to_data(debug_graph, "debug_data_output")
	debug_seq = dump_to_sequencegraphswithoutbonds("debug_dump_output")
	sequencegraphs_to_dump(debug_seq, "debug_seq_dump_output")
	dumpbonds_to_brokenbondssets("debug_dumpbd")
	renumber_id_data("debug_data_output")
	renumber_id_dump("debug_dump_output")

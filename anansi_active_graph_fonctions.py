#!/usr/bin/python -i

"""
---------------------------------------------------------
					Anansi 1.0 17/07
------------------------------------------------------------
Code for creation and caracterisation of elastomer networks
Coded by Morgane Mahaud for Mateis lab
Contact: morganemahaud@hotmail.com

Versions: python 2.7
------------------------------------------------------------
"""

__author__ = "Morgane Mahaud <morganemahaud@hotmail.com>"
__date__ = "31 July 2017"

__version__ = "$Revision: 1 $"

import random
import copy
import sys
sys.path.insert(1,'./..')
import networkx as nx
import numpy as np
from anansi_passive_graph_fonctions import 
from anansi_file_conversion import

	
#def selectionner_groupe_possible_prec
#autre solution: jouer avec attribut masque des noeuds
   
def compensate_spatial_drift(sequence):
	"""([nx.Graph]) -> [nx.Graph]
	Modify the positions of the units in a graph sequence so that any unwanted global linear movement disappear
	"""
	graph0 = sequence[0]
	nodes = graph0.nodes()
	pos0 = graph0.node[nodes[0]]["pos"]
	for graph in sequence[1:]:
		driftpos = graph.node[nodes[0]]["pos"] - pos0#the position without drift. To be more exact, the speed of unit 0 should be taken into account but is insignificant next to the drift
		for unit in nodes[1:]:
			graph.node[unit]["pos"] = graph.node[unit]["pos"] - driftpos
	return sequence

def count_units_by_attribut(graph, att = ()):#!!!!!
	"""(nx.Graph, tuple(str))-> dictionary	
	Return a dictonnary whose keys are attributes (ex: "type", "mol", "pos") and values are the number of nodes presenting this attribut
	"""
	dic = {att:None}
	for unit in graph:
		tmp = tuple([graph.node[unit][a] for a in att])
		if tmp in dic: dic[tmp] += 1
		else: dic[tmp] = 1
	return dic

def cut(graph, dc, zone = None, tobecut = None, sauvdata = None):
	"""(nx.Graph, float, [int], [(int,int)], bool, str) -> nx.Graph
	Cut bonds in a graph.
	
	Parameters:
	-------------
	graph: nx.Graph: graph to be modified
	dc: float: cut density, number of bonds cut divided by number of units in the box
	zone: [int] default = None: list of units from which bonds can be cut
	tobecut: [(int,int)] default = None: list of bonded that must be cut in priority
	sauvdata: string default = None: 
	
	Output:
	-------------
	res: nx.Graph: graph after cut
	"""
	print "Cut beginning..."
	compt = 0
	
	#selection of bonds to be cut, in the order given by tobecut (if not None) or among the zone (if not None). Else random among all bonds
	print "...selection of bonds to be cut"
	aim = int(dcl*len(graph))
	if tobecut:
		aimfrac = [int(i[0])*len(i[1]) for i in tobecut]
		aimfrac.append(aim - sum(aimfrac))
		assert aimfrac[-1] >= 0
		indaim = 0
		precs.append((float(aimfrac[-1]/len(graph)), graph.bonds()))
		print 'TO BE CUT:', [(i, len(j)) for i,j in precs]
		precurseurs = random.sample([i for i in precs[0][1]], aimfrac[0])
	elif zone:
		precurseurs = random.sample([i for i in zone], aim)
	else:
		precurseurs = random.sample([i for i in graph.edges()], aim)
	res = graph	
	#coupure des liens
	for a,b in precurseurs:	
		possible = True	
		#si possible, coupure du lien
		if possible:
			res.remove_edge(a, b)
			compt += 1
			res.node[a]['type'] = 4
			res.node[b]['type'] = 4
				
		#selection de nouveaux precurseurs si le nombre de liens aim n'est pas atteint
		if (a,b) is precurseurs[-1] and compt[0] < aim:
			if precs:
				if compt == sum(aimfrac[:indaim+1]):
					indaim += 1
					echant = aimfrac[indaim]
				else: echant == sum(aimfrac[:indaim+1]) - compt[0]
				precurseurs.extend(random.sample([n for n in precs[indaim] if not raph[n[0]][n[1]]['masque']], echant))
			if zone:
				precurseurs.extend(random.sample([n for n in zone if not raph[n[0]][n[1]]['masque']], aim-compt[0]))
			else:
				precurseurs.extend(random.sample([n for n in graph.edges() if not graph[n[0]][n[1]]['masque']], aim-compt[0]))
	print "...fin des coupures"
	res.graph['bonds'] = res.number_of_edges()
	if sauvdata: graph_to_data(res, sauvdata)
	return res
  
def more_bondtype(graph, indmax = None):
	if not indmax: indmax = graph.graph['bond types'] + 1
	bondcoeff = graph.graph['bond coeffs'][0][1:]
	for i in range(graph.graph['bond types']+1, indmax+1): 
		graph.graph['bond coeffs'].append(str(i)+bondcoeff)
	graph.graph['bond types'] = max(indmax, graph.graph['bond types'])

def particularites_to_type(graph, part, bonds = False, seuil = 3, bd2eqclk = True):
	"""particularites_to_type(graph, part, bonds = True) -> rien
	(nx.Graph, list ou string, boolean) -> None
	Change le type des atomes ou des liens ayant la/les particulatires part d'apres l'index
	Si bonds est vrai, le type des liens est change, sinon le type des atomes est change
	Si le lien a plusieurs particularites, la premiere de part sera retenue
	INDEX={"intra":4, "extra":5, "miniboucle":6, "moignon":7, "multilien":8, "doublet":9, "hexa":10, "moyeu":11}
	"""
	decompter_particularites(graph, seuil, bd2eqclk)
	if isinstance(part,str): part = [part]
	#creation des types de liens ou atomes supplementaires
	indmax = max([INDEX[p] for p in part])
	if bonds: more_bondtype(graph, indmax)
	else:
		#attention, probablement foireux
		paircoeff = graph.graph['pair coeffs'][0][3:]
		for i in range(graph.graph['atom types']+1, indmax+1):
			for j in range(1, i):
				graph.graph['pair coeffs'].append(str(j)+" "+str(i)+paircoeff)
		graph.graph['atom types'] = max(indmax, graph.graph['atom types'])
	for a,b in graph.edges_iter():
		if "part" in graph[a][b] and graph[a][b]['part']:
				p = [i for i in part if i in graph[a][b]['part']]
				if p: 
					if bonds: graph[a][b]['type'] = INDEX[p[0]]
					else:
						graph.node[a]['type'] = INDEX[p[0]]
						graph.node[b]['type'] = INDEX[p[0]]
						
def rename_id_by_molecule_order(graph):
	"""rename_id_by_molecule_order(graph) -> null
	Modified the units id in the graph so that for M chains of lenght Nj, the ist unit of the jst molecule have for id i+sum(k=0,j-1, Nk)
	/!\ If apply on elastomer
	"""
	listmolecules = get_sections(graph)
	M = len(listmolecules)
	sumN = 0
	dictmolecules = {graph.node[molecule[0]]['mol']: molecule for molecule in listmolecules}
	for j in range(1, M):
		N = len(dictmolecules[j])
		for i in range(0, N):
			
		sumN += N
			

def reticuler(graph, conds, dr, rmax = 1.5, seuil = 3, zone = None, precs = None, bondtype = 2, maxbd = 3, proximum = False, visuel = True, sauvdata = None, blabla=True):
	"""reticuler(graph, conds, dr, rmax, **kwargs) -> res
		(nx.Graph, [str], float, float, **kwargs) -> [nx.Graph]
	Cree des reticulations dans un graph de polymere ou elastomere
	Ajoute un type de lien designant les reticulations d'identifiant et coefficients reglables
	Passe les atomes lies par reticulation au type 1
	args:
		- graph: graph sur lequel travailler
		- conds: liste des particularites topologiques a retirer de l'elastomere. cf suite
		- dr: nombre de liens a cree/ nombre d'unites dans la boite
		- rmax: rayon de 'captation' des voisins
		- seuil = 3: int, nombre de liens en dessous duquel on considere une particularite topologique
		- zone = None: liste des identifiants d'unites pouvant reticuler entre elles, par defaut toute la boite
		- precs = None: [(float, [])], liste des selections strictes dans l'ordre de priorite
		- bondtype = None: int, type des liens crees. Par defaut, maximum des liens existant+1
		- maxbd = 3: int, nombre de liens maximal tolere par unite
		- proximum = False: bool, si True, le plus proche voisin dans la sphere de captation est selectionne, si False, l'un des voisins
		- visuel = True: bool, visualisation de la repartition des sections
		- sauvdata = None: str, creation d'un fichier data de lammps de nom sauvdata
	
	conds: 
	Conds se presente sous la forme ['a_b', 'c'...] ou a, b et c sont des conditions donnees si dessous.
	'c' donnera un reseau subissant la contrainte c, 'a_b' donnera un reseau subissant les contraintes a et b
	Le lien l ne sera pas cree si:
		- alea ou '': le precurseur ne trouve pas de cible (non compabilise)
		- doublet: il existe au moins 1 lien entre les 'seuil'voisins du precurseur et ceux de sa cible
		- extra: l relie 2 molecules differentes
		- hexa: il existe au moins 1 lien depuis les 'seuil'voisins du precurseur ou de sa cible vers l'exterieur
		- intra: l relie 2 unites de la meme molecule
		- miniboucle: l relie des unites reliees par seuil liens ou moins
		- moignon: l cree un bout de chaine de nombre de liens < seuil
		- moyeu: il existe au moins 1 miniboucle dans les 'seuil'voisins du precurseur ou de sa cible
		- multilien: il existe au moins 1 lien dans les 'seuil'voisins du precurseur ou de sa cible
	Le nombre prc*nbunites de liens a atteindre sera atteint par le premier reseau de la liste.
	Les reseaux suivant presenteront plus ou moins de liens que le premier si ses conditions sont differentes.
	Un decompte du nombre de liens evites pour chaque condition est donne. Ce decompte peut compter plusieurs
	fois le meme lien s'il obei a plusieurs conditions.
	
	mapinit:
	Conditions de choix des unites pouvant se lier. A ecrire
	
	"""
	if blabla: print "\nDebut de la reticulation..."
	if isinstance(conds, str): conds = [conds]
	if blabla: print "...restrictions appliquees:\n\t- "+"\n\t- ".join(conds)
	
	condsinst = list()
	comptsans = list()
	for sans in conds:
		if sans == 'alea': 
			condsinst.append(list())
			comptsans.append(dict())
		else: 
			condsinst.append([str_to_part(i) for i in sans.split('_')])
			comptsans = [{i: 0 for i in sans.split('_')}]
	nbg = len(conds)
	rgnbg = range(nbg)
	res = [graph]
	compt = [0 for i in rgnbg]
	dbimp = [0 for i in rgnbg]
	
	more_bondtype(graph)
	b = boxparam(graph)
	comparts = cadrillage(graph, b)
	
	#initialisation de la map des unites pouvant etre liees
	for at in graph:
		if not graph.node[at]['masque']: graph.node[at]['masque'] = graph.degree(at) >= maxbd

	#creation de copies de graphs pour les differentes conditions
	if blabla: print "...duplication du reseau"
	for i in range(nbg-1): res.append(graph.copy())
	
	#choix des precurseurs, dans l'ordre donne par precs, complete par aleatoire
	if blabla: print "...choix des precurseurs"
	if zone: cible = int(dr*len(zone))
	else: cible = int(dr*len(graph))
	if precs:
		ciblefrac = [int(i[0])*len(i[1]) for i in precs]
		ciblefrac.append(cible - sum(ciblefrac))
		assert ciblefrac[-1] >= 0
		indcible = 0
		precs.append((float(ciblefrac[-1]/len(graph)), graph.nodes()))
		print 'PRECS:', [(i, len(j)) for i,j in precs]
		precurseurs = random.sample([i for i in precs[0][1] if not graph.node[i]['masque']], ciblefrac[0])
	elif zone:
		precurseurs = random.sample([i for i in zone if not graph.node[i]['masque']], cible)
	else:
		precurseurs = random.sample([i for i in graph.nodes() if not graph.node[i]['masque']], cible)
		
		
	if blabla: print "...creation des liens"
	maxbd -= 2 #moins les deux voisins dans la molecule
	for tag in precurseurs:	
		if not graph.node[tag]['masque']:
			voisins = vois_zone(tag, rmax, graph, comparts, b, False, False, proximum)
			if not proximum and voisins: vois = random.choice(voisins)
			else: vois = voisins
			# creation ou non du lien dans les sorties avec leurs conditions propres
			if vois:
				for i in rgnbg:
					#verification de l'environnement du lien
					if conds[i]:
						possible = True
						envir = EnvirCrosslink(tag, vois, res[i], seuil)
						for instpart in condsinst[i]:
							if instpart.estPart(res[i], envir):
								if not possible: dbimp[i] += 1
								possible = False
								comptsans[i][instpart.nom()] += 1	
						if zone: possible = possible and vois in zone
					#si possible, creation du lien
					if possible:
						res[i].add_edge(tag, vois, type = bondtype)
						res[i].node[tag]['masque'] = True
						res[i].node[vois]['masque'] = True
						compt[i] += 1
						res[i].node[tag]['type'] = 1
						res[i].node[vois]['type'] = 1
				
			#selection de nouveaux precurseurs si le nombre de liens cible n'est pas atteint
			if tag is precurseurs[-1] and compt[0] < cible:
				if precs:
					if compt[0] == sum(ciblefrac[:indcible+1]):
						indcible += 1
						echant = ciblefrac[indcible]
					else: echant == sum(ciblefrac[:indcible+1]) - compt[0]
					precurseurs.extend(random.sample([n for n in precs[indcible] if not graph.node[n]['masque']], echant))
				if zone:
					precurseurs.extend(random.sample([n for n in zone if not graph.node[n]['masque']], cible-compt[0]))
				else:
					precurseurs.extend(random.sample([n for n in graph.nodes() if not graph.node[n]['masque']], cible-compt[0]))
					
	if blabla: print "...%i precurseurs supplementaires selectionnes" %(len(precurseurs)-cible)
	for i in rgnbg: 
		if blabla: print "...reseau %s: creation de %i liens" %(conds[i], compt[i])
		res[i].graph['bonds'] += compt[i]
		for item in comptsans[i].items(): print "\t - %i liens %s evites" %(item[1], item[0])
		if blabla: print "\t - %i surplus de conditions" %(dbimp[i])
	if blabla: print "...fin de la reticulation"
		
	if visuel:
		style = ['r-', 'b-', 'g-', 'k-', 'c-', 'm-', 'y-', 'r--', 'b--', 'g--', 'k--', 'c--', 'm--', 'y--']
		entree = list()
		leg = list()
		if blabla: print "\nNombres de pendants et moyennes des sections:"
		for i in rgnbg:
			rep = repartition_sections(res[i])
			rep[1] += nombre_singlets(res[i])-compt[i]
			nbpdt = len([0 for j in res[i].nodes_iter() if res[i].degree(j) == 1])
			msections = sum(rep)/len(rep)
			print "...%s: %i, %f" %(conds[i], nbpdt, msections)
			if conds[i] == 'alea': leg.append('alea')
			else: leg.append("exclu "+conds[i])
			entree.extend([range(len(rep)), rep, style[i]])
		f = plt.figure()
		ax = f.add_subplot(1,1,1)
		ax.plot(*entree)
		ax.legend(leg)
		ax.set_xlabel('segment repartition')
		ax.set_ylabel('nb of segments')
		plt.show()	
		
	if sauvdata:
		if "*" in sauvdata:
			split = sauvdata.split('*')
			for i in rgnbg: 
				if conds[i] == "": conds[i] = 'alea'
				graph_to_data(res[i], conds[i].join(split))
		elif nbg == 1: graph_to_data(res[i], sauvdata)
		else:
			for i in rgnbg: graph_to_data(res[i], sauvdata+'_'+conds[i])
	
	if len(conds) == 1: return res[0]
	return res

def retirer_zone(graph, zone):
	"""retirer_zone(graph, zone)-> None
	(nx.Graph, [int]) -> None
	Retire les atomes et liens du graph se trouvant dans la liste d'indentifiants donnee en entree
	"""
	graph.remove_nodes_from(zone)
	graph.graph['atoms'] = len(graph)
	graph.graph['bonds'] = nx.number_of_edges(graph)
	
def suivre_particularites_from_dumps(dump, dumpbd, part, nomsortie = None, seuil = 3, bd2eqclk = True, rompu = True):
	"""marquer_particularites_dump(dump, datainit, part, nomsortie = None, seuil = 3, bd2eqclk = True) -> None
		(str, str, [str], str, int, bool) -> None
	A partir d'un fichier dump de Lammps, marque les particularites demande en changeant le type d'atome d'apres l'index
	Ecrit un fichier nomsortie ou reecrit le fichier une fois les modifications faites
	Si bonds est vrai, le type des liens est change, sinon le type des atomes est change
	Si le lien a plusieurs particularites, la premiere de part sera retenue
	INDEX={"intra":4, "extra":5, "miniboucle":6, "moignon":7, "multilien":8, "doublet":9, "hexa":10, "moyeu":11}
	"""
	seq = dumps_to_seqbd(dump, dumpbd)
	for n in range(len(seq)):
		marquer_particularites(seq[n], seuil, bd2eqclk)
		print "--------------------------------------\nt = ", seq[n].graph['timestep']
		for p in part:
			for tag in seq[n].nodes(): 
				if not rompu or seq[n].node[tag]['type']!=INDEX['rompu']: seq[n].node[tag]['type'] = INDEX[p[0]]
	seq_to_dump(seq, nomsortie)

def suivre_rupture_from_data(dump, datainit, nomsortie = None, rmax = 1.5):
	"""
	"""
	initg = data_to_graph(datainit)
	edges = initg.edges()
	sqrmax = rmax**2
	seq = dump_to_seq(dump)
	print "Suivi de rupture..."
	comptot = 0
	for n in range(len(seq)):
		compt = 0
		boxhi = np.array([initg.graph['xhi'], initg.graph['yhi'], initg.graph['zhi']])
		boxlo = np.array([initg.graph['xlo'], initg.graph['ylo'], initg.graph['zlo']])
		box = boxhi-boxlo
		for a,b in edges:
			posb = compenser_periode_box(seq[n].node[a]['pos'], seq[n].node[b]['pos'], box)
			if sqdist(seq[n].node[a]['pos'], posb) > sqrmax: 
					seq[n].node[a]['type'] = INDEX['rompu']
					seq[n].node[b]['type'] = INDEX['rompu']
					compt += 1
		comptot += compt
		if compt: print "\tt= %s: %i liens rompus depuis le dernier pas, %i au total" %(seq[n].graph['timestep'], compt, comptot)
	print len([0 for i in seq[0].nodes() if seq[n].node[i]['type'] == 4])
	print "...fin du suivi de rupture."
	seq_to_dump(seq, nomsortie)
	
def suivre_rupture_from_dumps(dumpid, dumpbd, nomsortie = None):
	"""
	"""
	if isinstance(dumpid, list):
		seq, diff = [], []
		for dpid in dumpid: seq.extend(dump_to_seq(dpid))
		for dpbd in dumpbd: diff.extend(dumpbd_to_bkbd_sets(dpbd))
	else:
		seq = dump_to_seq(dumpid)
		diff = dumpbd_to_bkbd_sets(dumpbd)
	liste = list()
	comptot = 0
	for n in range(len(seq)):
		liste.extend(diff[n])
		for a, b in liste:
			seq[n].node[a]['type'] = INDEX['rompu']
			seq[n].node[b]['type'] = INDEX['rompu']
		compt = len(diff[n])
		comptot += compt
		if diff[n]: print "\tt= %s: %i liens rompus depuis le dernier pas, %i au total" %(seq[n].graph['timestep'], compt, comptot)
	seq_to_dump(seq, nomsortie)
	   
def suivre_zones_from_dump(dumpid, zones, nomsortie = None):
	"""suivre_zones_from_dump(dumpid, zonesid, nomsortie = None) -> None
		(str, list([int]), str) -> None
	"""
	print "Suivi de zones..."
	if not nomsortie: nomsortie = dumpid
	seq = dump_to_seq(dumpid)
	#creation du dictionnaire d'association id-mol
	dic = {}
	for z in range(len(zones)):
		dic.update([(i, z+2) for i in zones[z]])
	#remplacement des mols
	for g in seq:
		for at in g.nodes():
			if at in dic: g.node[at]['mol'] = dic[at]
			else: g.node[at]['mol'] = 1
	print "...fin du suivi de zones."
	seq_to_dump(seq, nomsortie)


if __name__ == '__main__':
	#routine pour verifier si toutes les fonctions repondent
	reticuler(debug_graph, ['alea'], 0.1, 1.5, visuel = False)
	particularites_to_type(debug_graph, ['intra'])

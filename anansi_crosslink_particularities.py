#!/usr/bin/python -

import sys


	
def str_to_particularity(s):
	"""(str) -> Particularity
	return a instance of the particularity of name s
	"""
	return getattr(sys.modules['anansi_crosslink_particularities'], s.capitalize())()

def get_kneighbors(graph, node, k, minusbond = False):
	"""(nx.Graph, Object, int, tuple) -> [int]
	return the list of id of nodes at a distance of k bonds or less from the node of id node in graph.
	node is included in the output list.
	If minusbond = (a, b), the bond (a, b) and its branchs are excluded from the output.
	"""
	kneighbors = [node]
	ind = 0
	for i in range(k): 
		lg = len(kneighbors)
		kneighbors += [v for n in kneighbors[ind:] for v in graph.neighbors(n) if v not in kneighbors and (n,v)!=minusbond and (v,n)!=minusbond]
		ind = len(kneighbors) - lg - 1
	return kneighbors
			
class BondSurrondings(object):
	"""Class describing the surroundings of a bond.

	ATTRIBUTS:
	a: int: id of the first unit belonging to the studied bond
	b: int: id of the second unit belonging to the studied bond
	akneighbors: [int]: list of kneighbors of a, excluding the ones from (a,b) branchs
	bkneighbors: [int]: list of kneighbors of b, excluding the ones from (a,b) branchs
	tria: [int]: list of kneighbors of a which belong to three bonds, excluding a and the ones from (a,b) branchs
	trib: [int]: list of kneighbors of b which belong to three bonds, excluding b and the ones from (a,b) branchs
	tri: [int]: list of kneighbors of a and b which belong to three bonds, excluding a and b
	"""

	def __init__(self, a, b, graph, k):
		"""__init__(self, a, b, graph, k) -> BondSurrondings
		create a BondSurrondings instance of bond (a, b) in the graph graph which include its kneighbors
		"""
		self.a = a
		self.b = b
		self.akneighbors = get_kneighbors(graph, a, k, (a,b))
		self.bkneighbors = get_kneighbors(graph, b, k, (a,b))
		self.tria = [i for i in self.akneighbors if graph.degree(i) == 3 and i!=a]
		self.trib = [i for i in self.bkneighbors if graph.degree(i) == 3 and i!=b]
		self.tri = self.tria+self.trib

class Particularity(object):
	"""Mother class describing a bond particularity
    """
	def __init__(self):
		pass

	def name(self):
		return self.__class__.__name__.lower()
		
	def index(self):
		return INDEX[self.name]
	
	def ispartic(self, graph, surround):
		pass
	
	def mark(self, graph, surround):
		"""mark two nodes as forming a particularity.
		The nodes get a attribut 'part' whose value is the name of the particularity 
		"""
		if self.ispartic(graph, surround): graph[surround.a][surround.b]['part'].append(self.name)
		
	
class Intra(Particularity):
	"""Class describing a bond linking two units of the same molecule/polymer chain
	"""
	def ispartic(self, graph, surround):
		return graph.node[surround.a]['mol'] == graph.node[surround.b]['mol']
	
class Extra(Particularity):
	"""Class describing a bond linking two units of different molecules/polymer chains
	"""
	def ispartic(self, graph, surround):
		return graph.node[surround.a]['mol'] != graph.node[surround.b]['mol']
	
class Tinycycle(Particularity):
	"""Class describing a bond forming a small cycle, i.e. linking two units of the same molecule/polymer chain distant of k bonds or less
	"""
	def ispartic(self, graph, surround):	
		return surround.b in surround.akneighbors
	
class Chainend(Particularity):
	"""Class describing a bond distant of k bonds or less from a chain segment end
	"""
	def ispartic(self, graph, surround):
		return any([graph.degree(i) == 1 for i in surround.akneighbors+surround.bkneighbors])
	
class Multibond(Particularity):
	"""Class describing a bond distant of k bonds or less from another crosslink
	"""
	def ispartic(self, graph, surround):	
		return any(surround.tri)
	
class Double(Particularity):
	"""Class describing a double bond, i.e. a bond distant of k bonds or less from another bond linking the same parts of chains
	"""
	def ispartic(self, graph, surround):	
		return any([j in graph.neighbors(i) for i in surround.tria for j in surround.trib])
	
class Hexa(Particularity):
	"""Class describing an hexavalent node, i.e. a bond distant of k bonds or less from another bond linking the same part of chain and another part of chain
	"""
	def ispartic(self, graph, surround):	
		return ((surround.tri and (not surround.tria or not surround.trib)) or\
			(not all([any([j in graph.neighbors(i) for j in surround.tria]) for i in surround.trib]) and\
			not all([any([j in graph.neighbors(i) for j in surround.trib]) for i in surround.tria])))


INDEX = {"broken":4, "intra":5, "extra":6, "tinycycle":7, "chainend":8, "multibond":9, "double":10, "hexa":11}
ALLPART = [i() for i in Particularity.__subclasses__()]
ALLPARTSTR = [i.name for i in Particularity.__subclasses__()]

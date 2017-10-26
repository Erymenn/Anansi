#!/usr/bin/python -i

import os, sys, copy
#sys.path.insert(1,'./..')
import networkx as nx
import numpy as np
from anansi_conversion_fichier import graph_to_dump#, debug_graph
from anansi_particularites_reticulation import *
from math import floor, ceil

#INDEX = {"rompu":4, "intra":5, "extra":6, "miniboucle":7, "moignon":8, "multilien":9, "doublet":10, "hexa":11, "moyeu":12}
sqdist = lambda a, b: sum((a-b)**2)

def bande(direc, lo, hi, graph, cadrillage, boxparam):
    """(direc, lo, hi, graph, cadrillage, boxparam) -> bande
    (int or str, float, float, nx.Graph, list(list(list([int]))), (np.array(float), np.array(float), np.array(float), np.array(int), np.array(int))) -> [int]
    Renvoie une liste de tous les atomes du graph compris entre les positions lo et hi dans la direction direc
    /!\ lo et hi doivent etre compris entre les valeurs inferieure et superieure de la boite dans la direction donnee
    """
    boxhi, boxlo, box, l, inf = boxparam
    dicdir = {'x':0, 'y':1, 'z':2}
    if isinstance(direc, str): direc = dicdir[direc]
    bande = list()
    zone = [range(l[0]),range(l[1]),range(l[2])]
    #ajout de tout les cubes de cadrillage qui sont dans la bande
    ilo = int(ceil(lo)) - inf[direc]
    ihi = int(floor(hi)) - inf[direc]
    zone[direc] = zone[direc][ilo:ihi]
    for kx in zone[0]:
        for ky in zone[1]:
            for kz in zone[2]: 
                bande.extend(cadrillage[kx][ky][kz])
    #ajout des atomes dans les cubes de cadrillage aux frontieres de la bande
    zone[direc] = (ilo-1, ihi+1)
    for kx in zone[0]:
        for ky in zone[1]:
            for kz in zone[2]:
                for at in cadrillage[kx][ky][kz]:
                    if graph.node[at]['pos'][direc] > lo and graph.node[at]['pos'][direc] < hi:
                        bande.append(at)
    return bande

def boxparam(graph):
    """boxparam(graph) -> boxhi, boxlo, box, l, inf
        (nx.Graph) -> np.array(float), np.array(float), np.array(float), np.array(int), np.array(int)
    De la boite d'un graph renvoie: les limites superieures, les limites inferieures, la taille, 
    la taille de l'enveloppe a l'entier superieur, les limites inferieures a l'entier superieur
    """
    boxhi = np.array([graph.graph['xhi'], graph.graph['yhi'], graph.graph['zhi']])
    boxlo = np.array([graph.graph['xlo'], graph.graph['ylo'], graph.graph['zlo']])
    box = boxhi-boxlo
    l = box.astype(int) + 2
    inf = (boxlo//1).astype(int)
    return boxhi, boxlo, box, l, inf

def cadrillage(graph, boxparam):
    """cadrillage(graph) -> cadrillage
        (nx.Graph) -> list(list(list([int])))
    Decoupe la boite en compartiments d'arrete 1
    Renvoie une matrice 3D des identifiants des atomes presents dans le compartiment
    Sert a reduire le nombre d'unites a tester pour trouver les plus proches voisins pour + d'efficacite
    Le graph doit comporter identifiants, positions et coordonnes des coins
    """
    boxhi, boxlo, box, l, inf = boxparam
    cadrillage = [[[[] for iz in xrange(l[2])] for iy in xrange(l[1])] for ix in xrange(l[0])]
    for at in graph:
        ix, iy, iz = (graph.node[at]['pos'].astype(int) - inf)%l
        cadrillage[ix][iy][iz].append(at)
    return cadrillage

def compenser_periode_box(posa, posb, box):
    """compenser_periode_box(posprec, poscour, box) -> poscour
        (np.array(float), np.array(float), np.array(float)) -> np.array(float)
    renvoie la position corrigee de poscour par rapport a posprec en tenant compte de la periodicite de la boite.
    """
    vec = posa - posb
    pbox = vec > 2*np.ones(3)
    mbox = vec < -2*np.ones(3)
    if any(pbox) or any(mbox):
        pbox = pbox.astype(int)
        mbox = mbox.astype(int)
        temp = tuple(posb)
        posb += (pbox-mbox)*box
    #vec2 = posa-posb
    #if any(vec2-vec): print sum(vec)**2, sum(vec2)**2
    return posb

def croix(graph, cadrillage, boxparam, f):
    boxhi, boxlo, box, l, inf = boxparam
    a = 1-(1-f)**(1./3)
    croix = set()
    centre = boxlo + box/2
    lo = centre - a*box/2
    hi = centre + a*box/2
    for i in xrange(3): 
        croix.update(bande(i, lo[i], hi[i], graph, cadrillage, boxparam))
    return list(croix)

def decompter_particularites(graph, seuil = 1, bd2eqclk = True):
    """decompter_particularites(graph, seuil = 1, bd2eqclk = True) -> None
        (nx.Graph, int, boolean) -> None
    Affiche le nombre de liens de chaque type de "defauts"
    Seuil est le nombre de liens en dessous duquel 2 reticulations sont comptees comme 1 noeud
    Si bd2eqclk = True, on indique que les liens de type 2 du graph d'entree correspondent aux reticulations (+ efficace)
    """
    print "Decompte des particularites..."
    cl = get_cross_links(graph, bd2eqclk)
    cpart = {k:0 for k in INDEX}
    for a,b in cl:
        graph[a][b]['part'] = list()
        envir = EnvirCrosslink(a, b, graph, seuil)
        for part in ALLPART:
            if part.estPart(graph, envir):
                cpart[part.nom()] += 1
                graph[a][b]['part'].append(part.nom())
    for it in cpart.iteritems():
        print "%s %i" %(it[0], it[1])

def ellipsoide(graph, cadrillage, boxparam, ax, ay, az, rVms = None, rVst = None, centre = None):
    """ellipsoide(graph, cadrillage, boxparam, x, y, z, rVms = None, centre = None) -> ellipsoide
        (nx.Graph, list(list(list([int]))), (np.array(float), np.array(float), np.array(float), np.array(int), np.array(int)), float, float, float, float, np.array(int)) -> [int]
    """
    boxhi, boxlo, box, l, inf = boxparam
    if centre == None: centre = boxlo + box/2
    else: centre = np.array(centre)
    if rVms or rVst:
        if rVms: rVst = 1./(1 + rVms)
        coeff = (3*box[0]*box[1]*box[2]*rVst/(4*np.pi*ax*ay*az))**(1./3)
        ax = coeff*ax
        ay = coeff*ay
        az = coeff*az
    sqa = np.array([ax**2, ay**2, az**2])
    ic = (centre.astype(int) - inf)%l
    etendu = int(max(ax,ay,az))+1
    infzone = ic-etendu
    supzone = ic+etendu
    zone = list()
    ellip = list()
    for kx in [i%l[0] for i in xrange(ic[0]-etendu, ic[0]+etendu+1)]:
        for ky in [i%l[1] for i in xrange(ic[1]-etendu, ic[1]+etendu+1)]:
            for kz in [i%l[2] for i in xrange(ic[2]-etendu, ic[2]+etendu+1)]:
                zone += cadrillage[kx][ky][kz]
    for at in zone:
        i = (graph.node[at]['pos'].astype(int) - inf)%l
        mbox = i>supzone
        pbox = i<infzone
        if any(pbox) or any(mbox):
            pbox = pbox.astype(int)
            mbox = mbox.astype(int)
            posat = graph.node[at]['pos'] + (pbox-mbox)*box
            if posat**2/sqa < 1: ellip.append(at)
        else:
            if graph.node[at]['pos']**2/sqa < 1: ellip.append(at)
    return ellip

def fraction_soluble(graph, tailles=True):
    gs = nx.connected_component_subgraphs(graph)
    l = [len(i) for i in gs]
    if tailles: print l
    return 1 - float(max(l))/len(graph)

def get_cross_links(graph, bd2eqclk = True):
    """get_cross_links(graph, bd2eqclk = True) -> liste des cross-links
        (nx.Graph, booleen) -> [tuple]
    Renvoie la liste des liens de reticulation d'un graph
    Si bd2eqclk est vrai, on indique au programme que les bond de type 2 sont la totalite des liens de reticulation. Sinon on ignore bond types.
    """
    print "Recherche des liens de reticulation..."
    #si les bonds 2 du fichier data existent et sont des liens de reticulation, on les lit en tant que tels
    #bd2eqclk = False
    if bd2eqclk:
        print "Lecture des liens de type 2..."
        res = [bd for bd in graph.edges_iter() if graph[bd[0]][bd[1]]['type'] == 2]
    else:
        #liens entre molecules
        res = [e for e in graph.edges_iter() if graph.node[e[0]]['mol']!=graph.node[e[1]]['mol']]
        print len(res)
        #simplification du graph
        sur = surgraph(graph, True)
        sur.remove_edges_from(res)
        #suppresion des liaisons bout de chaine-autre mol
        #(env 15 pr 500 chaines de 200 et 2000 liens, soit 1.5prc des bouts)
        for e in sur.edges_iter():
            if graph.node[e[0]]['mol'] != graph.node[e[1]]['mol']: sur.remove_edge(e[0],e[1])
        #separation des molecules
        mols = nx.connected_component_subgraphs(sur)
        compt = 0
        for mol in mols:
            nodes = mol.nodes()
            l = len(nodes)
            rl = range(l)
            #creation des matrices de travail
            mat = np.array(nx.adjacency_matrix(mol))
            temp = len(np.nonzero(mat)[0]) - 2*(l-1)
            compt += temp
            marq = np.zeros((l,l))
            indices = {nodes[i]: i for i in rl}
            #boucle avant gestion des indeterminations
            diff = True
            alire = set(rl)
            while diff:
                diff = False
                resolus = set()
                for n in alire:
                    nz = np.nonzero(mat[n])[0]
                    lnz = len(nz)
                    ones = [i for i in nz if mat[n,i] == 1]
                    lno = lnz-len(ones)
                    if lno > 2 or lnz > 3 or lno > lnz:
                        print "ERREUR ligne %i: %i chemins et %i chemins longs" %(n, lnz, lno)
                    #cas avec 2 chemins longs
                    elif lno == 2:
                        diff = True
                        resolus.add(n)
                        for v in ones:
                            marq[n,v] = -1
                            marq[v,n] = -1
                    #cas du bout connu
                    elif lnz == 1:
                        diff = True
                        resolus.add(n)
                        if mat[nz[0],n] == 1: marq[nz[0],n] = 1
                    #cas lno = 0 ou 1 et lnz = 2 ou 3
                    else:
                        mo = [i for i in nz if marq[n,i] == 1 or mat[n,i] > 1]
                        mm = [i for i in nz if marq[n,i] == -1]
                        nm = [i for i in nz if i not in mo and i not in mm]
                        lnm = len(nm)
                        if lnm == 0:
                            diff = True
                            resolus.add(n)
                        elif lnm == 1 and lnz == 3:
                            diff = True
                            resolus.add(n)
                            marq[n,nm[0]] = -1
                            marq[nm[0],n] = -1
                        #cas (lnm = 2 et lnz = 3) ou (lnm = 1 et lnz = 2) 
                        
                alire = alire.difference_update(resolus)
                if not alire: diff = False
            if alire: print "INDETERMINATION"
            #recuperation de tous les liens grace a marq
            for n in rl:
                for v in xrange(n,l):
                    if marq[n,v] == -1:
                        res.append((nodes[n], nodes[v]))
        print "...%i liens de reticulation listes dont %i intramoleculaires" %(len(res), compt)
    return res

def get_densite_actives_chimiques(graph, db=False):
    if not db: g = get_graph_sans_solubles(graph)
    else: g = graph
    retirer_pendants(g)
    nb = float(len(get_sections(g)))
    return nb/volume(graph)

def get_mean_active_functionnality(graph):
	"""(nx.Graph) -> float
	return the mean functionnality of the active network
	"""
	overgraph = surgraph(get_graph_actif(graph), hybride = False)
	return float(sum([overgraph.degree(i) for i in overgraph]))/len(overgraph)

def get_graph_actif(graph):
    """get_graph_actif(graph) -> graph
    get_graph_actif(nx.Graph) -> nx.Graph
    renvoie le graph initial prive de ses parties solubles et pendantes
    """
    g = get_graph_sans_solubles(graph)
    retirer_pendants(g)
    print "densite de chaines actives:", float(len(get_sections(g)))/volume(graph)
    return g
    
def get_graph_sans_solubles(graph):
    """get_graph_sans_solubles(graph) -> graph
    get_graph_sans_solubles(nx.Graph) -> nx.Graph
    renvoie le graph initial sans ses fractions solubles
    """
    gs = nx.connected_component_subgraphs(graph)
    l = [len(i) for i in gs]
    ind = l.index(max(l))
    return gs[ind]

def get_pendants(graph):
    """
    """
    g = graph.copy()
    res = []
    bouts = [i for i in g.nodes_iter() if g.degree(i) == 1]
    nds = set(i for i in g.nodes_iter() if g.degree(i) > 2)
    #print "avant: ", len(bouts)
    if bouts:
        if len(bouts)>1: fin = bouts[-1]
        else: fin = None
        for i in bouts:
            #print i
            if i in g:
                prec = None
                cour = i
                while cour not in nds:
                    if fin and cour == fin: break #cas ou il n'y a pas de fraction active
                    res.append(cour)
                    suiv = [j for j in graph.neighbors(cour) if j != prec][0]
                    #print "rem", cour, suiv
                    prec = cour
                    cour = suiv
        #print "apres: ", len([i for i in graph.nodes() if graph.degree(i) == 1])
        g.remove_nodes_from(res)
        res.extend(get_pendants(g))
    return res

def get_pendants2(graph):
	g=graph.copy()
	retirer_pendants(g)
	return set(graph.nodes()).difference(set(g.nodes()))

def get_sections(graph, **kwargs):
    """get_sections(graph) -> sections
        (nx.Graph) -> liste
    renvoie une liste des sections > 1, sous la forme (1e unite, derniere unite, longeur totale section)
    **kwargs peut etre un booleen de cle data ou n'importe quel attribut associe a un noeud du graph.
    Si la cle est data, renvoie une liste des sections sous comme tuple des identifiants des unites.
    Si la cle est un attribut, renvoie une liste des sections comme tuples des attributs des unites.
    
    ex: graph = 0-1-2-3-4-5-6-7                sections = [(0,2,2), (3,6,3), (8,10,2), (6,16,4)]
                    | |     |                
               8-9-10-11-12 13-14-15-16 
    
    si data, sections = [(0,1,2), (3,4,5,6), (8,9,10), (6,13,14,15,16)]
    """
    #print "Determination des sections..."
    #boucle de reconstruction des sections de longueurs 2 qui ne soit pas des points de reticulation
    impairs = [i for i in graph if graph.degree(i) == 1 or graph.degree(i) == 3]
    subg_sec2 = nx.subgraph(graph, impairs)
    noncrosslinks = [e for e in subg_sec2.edges_iter() if subg_sec2[e[0]][e[1]]['type']==1]
    
    duo = [i for i in graph if graph.degree(i) == 2]
    subg_sec = nx.subgraph(graph, duo)
    isoles = [i for i in subg_sec.nodes_iter() if subg_sec.degree(i) == 0]
    bouts = [i for i in subg_sec.nodes_iter() if subg_sec.degree(i) == 1]
    #boucle de reconstruction de sections longues. 
    #Comme on a que des portions de graph lineaires, les sections sont reconstitues de voisin en voisin
    #Un algorithme utilisant nx.has_path et nx.shortest_path a ete teste. 
    #Il est 100x plus lent pour des polymeres et pire encore pour des elastomeres
    vus = set()
    att = False
    data = 'data' in kwargs.keys()
    boxhi = np.array([graph.graph['xhi'], graph.graph['yhi'], graph.graph['zhi']])
    boxlo = np.array([graph.graph['xlo'], graph.graph['ylo'], graph.graph['zlo']])
    box = boxhi-boxlo
    if not data and kwargs.keys(): att = kwargs.keys()[0]
    #premiere partie de boucles pour les sections de longeur 3
    if data: 
        sections = noncrosslinks
        sections.extend([(graph.neighbors(i)[0], i, graph.neighbors(i)[1]) for i in isoles])
    elif att: 
        sections = [(graph.node[e[0]][att], graph.node[e[1]][att]) for e in noncrosslinks]
        sections.extend([(graph.node[graph.neighbors(i)[0]][att], graph.node[i][att], graph.node[graph.neighbors(i)[1]][att]) for i in isoles])
    else: 
        sections = [(e[0], e[1], 1) for e in noncrosslinks]
        sections.extend([tuple(graph.neighbors(i)+[2]) for i in isoles])
    for i in bouts:
        if i not in vus:
            prec = i
            cour = subg_sec.neighbors(i)[0]
            #ajout du noeud a 1 ou 3 degree du debut de section
            a = [k for k in graph.neighbors_iter(i) if k != cour][0]
            if data: sec = [a, prec, cour]
            elif att: 
                sec = [graph.node[a][att]]
                sec.append(compenser_periode_box(sec[0], graph.node[prec][att], box))
                sec.append(compenser_periode_box(sec[1], graph.node[cour][att], box))
            else: lg = 1
            while subg_sec.degree(cour)!=1:
                tmp = cour
                cour = [k for k in subg_sec.neighbors_iter(cour) if k != prec][0]
                prec = tmp
                if data: sec.append(cour)
                elif att: 
                    if att == 'pos': 
                        sec.append(compenser_periode_box(sec[-1], graph.node[cour][att], box))
                    else: sec.append(graph.node[cour][att])
                else: lg += 1
            vus.add(cour)
            b = [k for k in graph.neighbors_iter(cour) if k != prec][0]
            if data: sections.append(sec+[b])
            elif att: 
                sections.append(sec+[compenser_periode_box(sec[-1], graph.node[b][att], box)])
            else: sections.append((a,b,lg+2))
    return sections

def get_solubles(graph):
    """get_solubles(graph) -> liste des sous-reseaux solubles
    get_solubles(nx.Graph) -> [[int]]
    renvoie une liste des liste d'atomes appartenant aux sous graphs non relies au reseau principal
    """
    gs = nx.connected_component_subgraphs(graph)
    l = [len(i) for i in gs]
    ind = l.index(max(l))
    gs.pop(ind)
    return [g.nodes() for g in gs]

def get_unfolded_box_sections(graph):
    """get_unfolded_box_sections(graph)-> sections
        (nx.Graph)-> [[np.array(3)]]
    Retour toutes les sections contenues dans une boite et etendues hors de la boite
    Si une section est coupee par un bord de boite elle est doublee et chacune prend les positions de chaque cote de la boite
    """
    boxhi, boxlo, box, l, inf = boxparam(graph)
    centralsecs = get_sections(graph, pos = True)
    centralbouts = [sec[i] for i in [0,-1] for sec in centralsecs]
    plussecs = []
    dirdeps = [np.array([x, y, z]) for x in (-1, 0, 1) for y in (-1, 0, 1) for z in (-1, 0, 1) if x or y or z]
    for d in dirdeps:
        for sec in centralsecs:
            incentral = False
            dsec = [pos + d*box for pos in sec]
            incentral = any([all(pos<boxhi) and all(pos>boxlo) for pos in dsec])
            if incentral and any([not any(dsec[i]-b) for i in [0, -1] for b in centralbouts]):
            #equiv a dsec[0] not in centralbouts and dsec[-1] not in centralbouts, puisque array ne peut etre lu comme booleen
                plussec.append(dsec)
    return centralsecs + plussecs        

def marquer_particularites(graph, seuil = 3, bd2eqclk = True):
    """
    """
    cl = get_cross_links(graph, bd2eqclk)
    for a,b in cl:
        graph[a][b]['part'] = list()
        envir = EnvirCrosslink(a, b, graph, seuil)
        for part in ALLPART:
            part.marquer(graph, envir)

def marquer_pendants_solubles(graph, dump = None, molactif = 0, db=False):
    """marquer_pendants_solubles(graph, dump = None, molactif = 0) -> None
    marquer_pendants_solubles(nx.Graph, str, int)-> None
    modifie l'identifiant de la molecule d'origine des segments de chaines du reseau pour la remplacer par molactif pour les chaines actives, molactif+1 pour les chaines pendantes et molactif+1+i pour les i composantes de la fraction soluble
    """
    size = graph.number_of_nodes()
    sol = get_solubles(graph)
    fs = float(sum([len(i) for i in sol]))/size
    print "fraction soluble:", fs
    if not db: pen = get_pendants(get_graph_sans_solubles(graph))
    else: pen = get_pendants(graph)
    fp = float(len(pen))/size
    print "fraction de chaines pendantes:", fp
    for n in graph.nodes_iter(): graph.node[n]['mol'] = molactif
    for n in pen: graph.node[n]['mol'] = molactif + 1
    i = 1
    if not db:
		for s in sol:
			i += 1
			for n in s: graph.node[n]['mol'] = molactif + i
    if dump:
        graph_to_dump(graph, dump)
    return fs, fp

def nombre_singlets(graph):
    """
    """
    trio = [i for i in graph if graph.degree(i) == 3]
    solo = [i for i in graph if graph.degree(i) == 1]
    clusters = nx.subgraph(graph, trio+solo)
    return nx.number_of_edges(clusters)
        
def order_atoms(fich):
	g = data_to_graph(fich)
	secs = get

def reconstituer_chemin_primitif(pos_sec):
	res = list()
	for sec in pos_sec:
		newsec = sec[:1]
		reste = set(sec[1:])
		for at in range(1, len(sec)):
			minsq = 4
			minpos = None
			for pos in reste:
				sqd = sqdist(newsec[-1], pos)
				if sqd < minsq:
					minsq = sqd
					minpos = pos
			newsec.append(minpos)
			reste.remove(minpos)
		res.append(newsec)
	return res
       
def repartition_sections(graph):
    """repartition_sections(graph) -> rep
        (nx.Graph) -> list
    Renvoie une liste ou rep[i] = nombre de sections de longueur i 
    """
    rep = []
    for s in get_sections(graph):
        if s[2] > len(rep)-1: rep.extend([0 for i in xrange(len(rep), s[2]+1)])
        rep[s[2]] += 1
    return rep

def repartition_sections_topologie(graph, nom, ttsol = False):
    """repartition_sections_topologie(graph) -> rep
        (nx.Graph) -> [list()]
    Renvoie un fichier des nombres de sections de longueur i pour l'ensemble, les parties actives, pendantes et solubles.
    Les molecules du fichier d'entree doivent correspondre aux-dites parties
    """
    rep = [[], [], [], []]
    for s in get_sections(graph, mol=True):
        mols = max(s[0], s[-1])+1
        lens = len(s)
        if ttsol:
            if mols > len(rep)-1: 
                for i in xrange(len(rep), mols+1): rep.append([])
        else:
            if mols > 2: mols = 3
        if lens > len(rep[0])-1: rep[0].extend([0 for i in xrange(len(rep[0]), lens+1)])
        if lens > len(rep[mols])-1: rep[mols].extend([0 for i in xrange(len(rep[mols]), lens+1)])
        rep[0][lens] += 1
        rep[mols][lens] += 1
    nbmol = len(rep)
    if not nom in os.listdir(os.getcwd()): open(nom, 'w').close()
    with open(nom,'w') as f:
        if ttsol: f.write("nb_d'unites total actif pendant soluble"+" soluble".join([str(i) for i in xrange(2, nbmol)])+'\n')
        else: f.write("nb_d'unites total actif pendant soluble\n")
        for j in xrange(1, max([len(i) for i in rep])):
            f.write(str(j))
            for i in xrange(nbmol):
                if j<len(rep[i]): f.write(" "+str(rep[i][j]))
                else: f.write(" 0")
            f.write('\n')
    return rep

def retirer_pendants(graph):
    """
    """
    for i in (j for j in graph.nodes_iter() if graph.degree(j) == 1):#for every chain end
           prec = None
           cour = i
           while graph.degree(cour) == 1:
               suiv = [j for j in graph.neighbors_iter(cour) if j != prec][0]
               graph.remove_node(cour)
               prec = cour
               cour = suiv

def sphere(graph, cadrillage, boxparam, f = None, centre = None, rayon = None, rVmh = None):
    """sphere(rVms, graph, cadrillage, boxparam, centre=None, rayon=None) -> sphere
        (float, nx.Graph, list(list(list([int]))), (np.array(float), np.array(float), np.array(float), np.array(int), np.array(int)), [int], float) -> [int]
    Attention!!! 
    si on utilise le rayon et centre pour definir la/les spheres, les intersections et auto-intersections ne sont pas prises en compte dans le calcul du volume
    si on utilise la fraction, l'auto-intersection due a la periodicite n'est prise en compte dans le volume que jusqu'a r<sqrt(2)/2 (soit f~0.965)
    """
    boxhi, boxlo, box, l, inf = boxparam
    if centre == None: centre = boxlo + box/2
    else: centre = np.array(centre)
    if rVmh: f = 1./(1 + rVmh)
    if not rayon: 
		r = (3*f/(4*np.pi))**(1./3)
		#if r > 0.5:#il faut retirer le volume des calottes qui depassent de la boite
			##f vaut alors pi/3*(-8*r**3+9*r**2-3./4)
			#coeffs = (8, -9, 0, 3./4+3./np.pi*f)
			#rs = [float(np.real(i)) for i in np.roots(coeffs) if np.isreal(i)]
			#rs2 =[i for i in rs if i>0]
			#r = min(rs2)
		rayon = r*(box[0]*box[1]*box[2])**(1./3)
    ic = (centre.astype(int) - inf)%l
    etendu = int(rayon)+1
    infzone = ic-etendu
    supzone = ic+etendu
    sqr = rayon**2
    zone = list()
    sphere = list()
    for kx in [i%l[0] for i in xrange(ic[0]-etendu, ic[0]+etendu+1)]:
        for ky in [i%l[1] for i in xrange(ic[1]-etendu, ic[1]+etendu+1)]:
            for kz in [i%l[2] for i in xrange(ic[2]-etendu, ic[2]+etendu+1)]:
                zone += cadrillage[kx][ky][kz]
    for at in zone:
        i = (graph.node[at]['pos'].astype(int) - inf)%l
        mbox = i>supzone
        pbox = i<infzone
        if any(pbox) or any(mbox):
            pbox = pbox.astype(int)
            mbox = mbox.astype(int)
            posat = graph.node[at]['pos'] + (pbox-mbox)*box
            if sqdist(centre, posat) < sqr: sphere.append(at)
        else:
            if sqdist(centre, graph.node[at]['pos']) < sqr: sphere.append(at)
    return sphere

def sphere_no_overlap(graph, cadrillage, boxparam, f = None, centre = None, rayon = None, rVmh = None):
    """sphere(rVms, graph, cadrillage, boxparam, centre=None, rayon=None) -> sphere
        (float, nx.Graph, list(list(list([int]))), (np.array(float), np.array(float), np.array(float), np.array(int), np.array(int)), [int], float) -> [int]
    Attention!!! 
    si on utilise le rayon et centre pour definir la/les spheres, les intersections et auto-intersections ne sont pas prises en compte dans le calcul du volume
    si on utilise la fraction, l'auto-intersection due a la periodicite n'est prise en compte dans le volume que jusqu'a r<sqrt(2)/2 (soit f~0.965)
    """
    boxhi, boxlo, box, l, inf = boxparam
    if centre == None: centre = boxlo + box/2
    else: centre = np.array(centre)
    if rVmh: f = 1./(1 + rVmh)
    if not rayon: 
		r = (3*f/(4*np.pi))**(1./3)
		if r > 0.5:#il faut retirer le volume des calottes qui depassent de la boite
			#f vaut alors pi/3*(-8*r**3+9*r**2-3./4)
			coeffs = (8, -9, 0, 3./4+3./np.pi*f)
			rs = (float(np.real(i)) for i in np.roots(coeffs) if np.isreal(i))
			rs2 =[i for i in rs if i>0]
			r = min(rs2)
		rayon = r*(box[0]*box[1]*box[2])**(1./3)
    ic = (centre.astype(int) - inf)%l
    etendu = int(rayon)+1
    infzone = ic-etendu
    supzone = ic+etendu
    sqr = rayon**2
    zone = list()
    sphere = set()
    for kx in (i%l[0] for i in xrange(ic[0]-etendu, ic[0]+etendu+1)):
        for ky in (i%l[1] for i in xrange(ic[1]-etendu, ic[1]+etendu+1)):
            for kz in (i%l[2] for i in xrange(ic[2]-etendu, ic[2]+etendu+1)):
                zone += cadrillage[kx][ky][kz]
    for at in zone:
        i = (graph.node[at]['pos'].astype(int) - inf)%l
        mbox = i>supzone
        pbox = i<infzone
        if any(pbox) or any(mbox):
            pbox = pbox.astype(int)
            mbox = mbox.astype(int)
            posat = graph.node[at]['pos'] + (pbox-mbox)*box
            if sqdist(centre, posat) < sqr: sphere.add(at)
        else:
            if sqdist(centre, graph.node[at]['pos']) < sqr: sphere.add(at)
    return list(sphere)

def sq_deplacement_moyen(graph0, graphi):
    """sq_deplacement_moyen(graph0, graphi)->distance**2
    sq_deplacement_moyen(nx.Graph, nx.Graph)->float
    renvoie le deplacement moyen des unites entre le graph i et le graph 0
    """
    res = [sqdist(graph0.node[n]['pos'], graphi.node[n]['pos']) for n in graph0.nodes_iter()]
    return sum(res)/len(res)

def surgraph(graph, hybride = False, seuil = None):
    """surgraph(graph, hybride = False, get_sections = False) -> surgraph
        (nx.Graph, bool, bool) -> nx.Graph ou nx.MultiGraph
    Si hybride == True, renvoie un graph dont les noeuds sont les noeuds a 1 ou 3 liaisons de graph (data inclu)
        et dont les liens sont les chemins entre ces noeuds, ponderes de la longueur du chemin
    Si hybride == False, renvoie un graph dont les noeuds sont les clusters de graph (data non inclu)
        et dont les liens sont les chemins entre ces clusters, ponderes de la longueur du chemin
        
    ex: graph = 0-1-2-3-4-5-6-7            hybride = 0-2-3-6-7        no hybribe = 0-(2,3,10,11,12)-(6,7)
                    | |     |                         / /   \                            |           |
               8-9-10-11-12 13-14-15-16           8-10-11-12 16                          8           16
               
        seuil=2 0-2-3  6-7
                  | |
               8-10-11-12
    """
    print "Creation du surgraph..."
    sections = get_sections(graph)    
    trio = [i for i in graph if graph.degree(i) == 3]
    solo = [i for i in graph if graph.degree(i) == 1]
    clusters = nx.subgraph(graph, trio+solo)
    if hybride:
        for a,b in clusters.edges_iter(): clusters[a][b]['weight'] = 1
        if seuil > 1: clusters.add_weighted_edges_from([s for s in sections if s[2] <= seuil])
        elif not seuil: clusters.add_weighted_edges_from(sections)
        return clusters
    else:
        dic = {}
        for n in clusters.nodes_iter():
            m = nx.node_connected_component(clusters, n)
            if len(m) == 1: m = n
            else:
                m.sort()
                m = tuple(m)
            dic[n] = m
        res = nx.MultiGraph()
        res.add_weighted_edges_from([(dic[s[0]], dic[s[1]], s[2]) for s in sections])
        return res
 
def unwrap(graph):
    """
    """
    vus = set()
    duo = [i for i in graph if graph.degree(i) == 2]
    subg_sec = nx.subgraph(graph, duo)
    boxhi = np.array([graph.graph['xhi'], graph.graph['yhi'], graph.graph['zhi']])
    boxlo = np.array([graph.graph['xlo'], graph.graph['ylo'], graph.graph['zlo']])
    box = boxhi-boxlo
    for i in (j for j in subg_sec.nodes_iter() if subg_sec.degree(j) == 1):#for every chain end
        if i not in vus:
            prec = i
            cour = subg_sec.neighbors(i)[0]
            preprec = [k for k in graph.neighbors_iter(i) if k != cour][0]
            temp = compenser_periode_box(graph.node[preprec]['pos'], graph.node[prec]['pos'], box)
            graph.node[prec]['pos'] = temp
            temp = compenser_periode_box(graph.node[prec]['pos'], graph.node[cour]['pos'], box)
            graph.node[cour]['pos'] = temp
            while subg_sec.degree(cour)!=1:
                tmp = cour
                cour = [k for k in subg_sec.neighbors_iter(cour) if k != prec][0]
                prec = tmp
                temp = compenser_periode_box(graph.node[prec]['pos'], graph.node[cour]['pos'], box)
                graph.node[cour]['pos'] = temp
            vus.add(cour)
            suiv = [k for k in graph.neighbors_iter(cour) if k != prec][0]
            temp = compenser_periode_box(graph.node[cour]['pos'], graph.node[suiv]['pos'], box)
            graph.node[suiv]['pos'] = temp

def vois_zone2(tag, graph, rmax, cadrillage, boxparam, masquevois = True, bdvois = True, proximum = False):
    """~~~
    """
    boxhi, boxlo, box, l, inf = boxparam
    zone = sphere(graph.node[tag], rmax, cadrillage, boxparam)
    for at in zone:
        visible = at != tag
        visible = visible and (bdvois or at not in graph[tag])
        visible = visible and (masquevois or not graph.node[at]['masque'])
        if visible:
            dist = np.linalg.norm(graph.node[at]['pos'] - graph.node[tag]['pos'])
            if dist < rmax:
                if proximum:
                    if dist < tmpdist: 
                        vois = at
                        tmpdist = dist
                else: voisins.append(at)
    if proximum: return vois
    return voisins

def vois_zone(tag, rmax, graph, cadrillage, boxparam, masquevois = True, bdvois = True, proximum = False):
    """vois_zone(tag, graph, cadrillage, rmax, masquevois = True, bdvois = True, proximum = False) -> voisin(s)
        (int, nx.Graph, [[[[int]]]], float, bool, bool, bool) -> [int] ou int
    Renvoie les voinsins contenus dans une sphere de centre tag et de rayon rmax, ou le voisin le plus proche
    Si masquevois = True, les voisins masques sont inclus
    Si bdvois = True, les voisins lies a tag sont inclus
    Si proximum = True, seul le voisin le plus proche est renvoye
    """
    #delimitation de la zone de recherche
    boxhi, boxlo, box, l, inf = boxparam
    etendu = int(rmax)+1
    ix, iy, iz = (graph.node[tag]['pos'].astype(int) - inf)%l
    zone = list()
    for kx in (i%l[0] for i in xrange(ix-etendu, ix+etendu+1)):
        for ky in (i%l[1] for i in xrange(iy-etendu, iy+etendu+1)):
            for kz in (i%l[2] for i in xrange(iz-etendu, iz+etendu+1)):
                zone += cadrillage[kx][ky][kz]
    #recherche des ou du plus proche(s) voisin(s)
    #y compris ceux deja lies a l'unite si bdvois et/ou ceux masques si masquevois 
    voisins = []
    vois = None
    tmpdist = rmax
    for at in zone:
        visible = at != tag
        visible = visible and (bdvois or at not in graph[tag])
        visible = visible and (masquevois or not graph.node[at]['masque'])
        if visible:
            dist = np.linalg.norm(graph.node[at]['pos'] - graph.node[tag]['pos'])
            if dist < rmax:
                if proximum:
                    if dist < tmpdist: 
                        vois = at
                        tmpdist = dist
                else: voisins.append(at)
    if proximum: return vois
    return voisins

def volume(graph):
    boxhi = np.array([graph.graph['xhi'], graph.graph['yhi'], graph.graph['zhi']])
    boxlo = np.array([graph.graph['xlo'], graph.graph['ylo'], graph.graph['zlo']])
    box = boxhi-boxlo
    return box[0]*box[1]*box[2]

#debug_pos_secs = get_sections(debug_graph, pos = True)    

if __name__ == '__main__':    
    #routine pour verifier si toutes les fonctions repondent
    #debug_graph = 0-1-2-3-4-5-6-7    
    #                      | |     |
    #              8-9-10-11-12 |  
    #                           13-14-15-16    
    #kvois(debug_graph, 1, 1)
    #nombre_singlets(debug_graph)
    #sq_deplacement_moyen(debug_graph, debug_graph)
    #compenser_periode_box(np.zeros(3), np.ones(3), 2*np.ones(3))
    #get_sections(debug_graph)
    #repartition_sections(debug_graph)
    #surgraph(debug_graph)
    #get_cross_links(debug_graph)
    #decompter_particularites(debug_graph)
    debug_cadrillage = cadrillage(debug_graph)
    vois_zone(3, debug_graph, debug_cadrillage, 2)

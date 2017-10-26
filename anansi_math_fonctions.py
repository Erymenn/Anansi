#!/usr/bin/python -i

import os
import numpy as np
import matplotlib.pyplot as plt
from anansi_passive_graph_fonctions import get_segments, boxparam, compenser_periode_box#, debug_pos_secs

sqdist = lambda a, b: sum((a-b)**2)

def mean_segments_path_lenght(pos_segments):
	"""[[np.array]] -> float
	return the average path lenght of the chain segments of a graph
	use pos_section, the list of the chain segments represented by the list of units positions (from get_segments)
	"""
	l = 0
	nb = 0
	for seg in pos_segments:
		prec = seg[0]
		for pos in seg[1:]:
			l += sum((pos - prec)**2)
			c += 1
			prec = pos
	return np.sqrt(l)/c

def barycentres(pos_segments):
	"""renvoie un tableau des coordonnees du barycentre de la molecule
	"""
	return [sum(seg)/len(seg) for seg in pos_segments]

def sq_dist_extrem_seg(pos_segments, liste = False):
	"""renvoie la distance quadratique moyenne entre extremites de segments
	"""
	if liste: return [sum((seg[-1]-seg[0])**2) for seg in pos_segments]
	return sum([sum((seg[-1]-seg[0])**2) for seg in pos_segments])/len(pos_segments)

def dist_extrem_seg(pos_segments, liste = False):
	if liste: return [np.sqrt(sum((seg[-1]-seg[0])**2)) for seg in pos_segments]
	return sum([np.sqrt(sum((seg[-1]-seg[0])**2)) for seg in pos_segments])/len(pos_segments)
	
def lpp(pos_segments, echant=1, liste = False):
	"""renvoie la longueur de contour moyenne des chaines
	"""
	res = []
	for seg in pos_segments:
		l=0
		prec = seg[0]
		i=0
		for pos in seg[1:]:
			if i == echant or not all(pos-seg[-1]):
				l += np.sqrt(sum((pos - prec)**2))
				prec = pos
				i=0
			i += 1
		res.append(l)
		#print res
	if liste: return res
	return sum(res)/len(res)

def moy_n(segments):
	"""Renvoie le nombre d'unites moyen des segments donnees en entree
	"""
	return sum([len(seg) for seg in segments])/len(segments)
	
#def l_persistance(graph):
	#s = get_segments(g, pos = True)
	#l = sum([for e in graph.edges_iter()])/graph.number_of_edges()
	#return "toto"
	
#def fractions_chaines_par_formules(graph, Ne=35):
	#"""formules de Charlesby_pinner, p97 these andre, polydisperse
	#return factives, fpend, fsol, facchim, facench
	#"""
	#nb = len(graph)
	#fs = fraction_soluble(graph)
	#sqfs = np.sqrt(fs)
	#d = 1#densite polym
	#q = 0.02#proba de retic
	#fc = q*d/2.*(6*sqfs*(1-sqfs)**3+2*(1-sqfs)**4)
	#fe = d/(2*Ne)*(2*(1-sqfs)**4))
	#return fe+fc, 1-fe-fc-fs, fs, fc, fe
	

#def fraction_actives_chimiques(graph)
	
def sq_rayon_giration(pos_segments, liste = False):
	g = barycentres(pos_segments)
	l = len(pos_segments)
	rgs = []
	for i in range(l):
		rg = [sum((pos-g[i])**2)  for pos in pos_segments[i]]
		rgs.append(sum(rg)/len(pos_segments[i]))
	if liste: return rgs
	return sum(rgs)/l

def sq_rayon_giration2(pos_segments, liste = False):
	l = len(pos_segments)
	rgs = []
	for i in range(l):
		rg = 0
		for j in pos_segments[i]:
			for k in pos_segments[i]:
				rg += (j-k)**2
		rgs.append(0.5*rg/(len(pos_segments[i])+1)**2)
	if liste: return rgs
	return sum(rgs)/l
	
def msid(pos_segments):
	summsid = []
	compt = []
	for seg in pos_segments:#iteration sur les molecules
		l = len(seg)
		lcompt = len(compt)
		summsid += [0 for i in xrange(lcompt, l)]
		compt += [0 for i in xrange(lcompt, l)]
		for n in xrange(1, l):
			temp = 0
			for i in xrange(l-n): #iteration dans la molecule
				plus = sum((seg[i+n]-seg[i])**2)#somme des coordonnees=distance au carre
				temp += plus
			summsid[n] += temp/(l-n)
			compt[n] += 1
	return [0]+[summsid[i]/compt[i] for i in xrange(1,len(compt))]

#def msid_and_err(pos_segments):
	#m = msid(pos_segments)
	#err = [0 for i in range(len(m))]
	#compt = [0 for i in range(len(m))]
	#for seg in pos_segments:#iteration sur les molecules
		#l = len(seg)
		#for n in range(1, l):
			#temp = 0
			#for i in range(l-n): #iteration dans la molecule
				#plus = (sum((seg[i+n]-seg[i])**2) - m[n])**2#somme des coordonnees=distance au carre - la moyenne pour n
				#temp += plus
			#err[n] += temp/(l-n)
			#compt[n] += 1
	#return m, [0]+[np.sqrt(err[i]/compt[i]) for i in range(1,len(m))]
	
def repartition_sid(pos_segments):
	nom = 'repartition.dat'
	res = {}
	for seg in pos_segments:
		l = len(seg)
		for n in xrange(1,l):
			for i in xrange(l-n):
				x = round(sum((seg[i+n]-seg[i])**2), 2)
				if x in res: res[x][n] += 1
				else: 
					res[x] = [0 for j in xrange(l)]
					res[x][n] += 1
	keys = res.keys()
	keys.sort()
	print "Debut de l'ecriture du fichier %s..." %(nom)
	if not nom in os.listdir(os.getcwd()): open(nom, 'w').close()
	with open(nom,'w') as f:
		f.write('x '+' '.join([str(i) for i in xrange(1, l)])+'\n')
		for r in keys:
			f.write(str(r)+' '+' '.join([str(res[r][i]) for i in xrange(1, l)])+'\n')
			
def repartition_err(pos_segments):
	m = msid(pos_segments)
	nom = 'repartition_err.dat'
	res = {}
	for seg in pos_segments:
		l = len(seg)
		for n in xrange(1,l):
			for i in xrange(l-n):
				x = round(sum((seg[i+n]-seg[i])**2) - m[n], 2)
				if x in res: res[x][n] += 1
				else: 
					res[x] = [0 for j in xrange(l)]
					res[x][n] += 1
	keys = res.keys()
	keys.sort()
	print "Debut de l'ecriture du fichier %s..." %(nom)
	if not nom in os.listdir(os.getcwd()): open(nom, 'w').close()
	with open(nom,'w') as f:
		f.write('x '+' '.join([str(i) for i in xrange(1, l)])+'\n')
		for r in keys:
			f.write(str(r)+' '+' '.join([str(res[r][i]) for i in xrange(1, l)])+'\n')
				
def msid_and_err(pos_segments):
	m = msid(pos_segments)
	err = [0 for i in xrange(len(m))]
	N = len(m)
	M = len(pos_segments)
	for n in xrange(1,N):
		c = 0
		for seg in pos_segments:
			c += 1
			for i in xrange(N-n):
				err[n] += (sum((seg[i+n]-seg[i])**2) - m[n])**2
		#print err[n], m[n], N-n, M,  (err[n]/((N-n)*M) - m[n]**2)/((N-n)*M)
		err[n] = np.sqrt(err[n])/((N-n)*M-1)
		#err[n] = np.sqrt((err[n]/((N-n)*M) - m[n]**2)/((N-n)*M))
	return m, err

#def msid_and_err(pos_segments):
	#m = msid(pos_segments)
	#err = [0 for i in range(len(m))]
	#N = len(m)
	#M = len(pos_segments)
	#for n in range(1,N):
		#c = 0
		#for seg in pos_segments:
			#c += 1
			#for i in range(N-n):
				#err[n] += (sum((seg[i+n]-seg[i])**2) - m[n])**2
				#if c<10: print sum((seg[i+n]-seg[i])**2) - m[n]
			#if c< 10: print n, c, err[n], np.sqrt(err[n]/((N-n)*c))
		#err[n] = np.sqrt(err[n]/((N-n)*M))
	#return m, err

#def msid_and_err(pos_segments):
	#m = msid(pos_segments)
	#err = [0 for i in range(len(m))]
	#compt = [0 for i in range(len(m))]
	#permol = [0 for i in range(len(m))]
	#x=0
	#for seg in pos_segments:#iteration sur les molecules
		#l = len(seg)
		#for n in range(1, l):
			#temp = 0
			#for i in range(l-n): #iteration dans la molecule
				#plus = sum((seg[i+n]-seg[i])**2)#somme des coordonnees=distance au carre
				#permol += plus
			#permol = permol/(l-n)
			#err[n] += (permol[n] - m[n])**2
			#compt[n] += 1
		#x+=1
		#if x<10:print err
	#print m, err, compt
	#return m, [0]+[np.sqrt(err[i]/compt[i]) for i in range(1,len(m))]

def see_msid(pos_segments, error = True):
	if error:
		m, err = msid_and_err(pos_segments)
		m2 = [m[i]/(m[1]*i) for i in xrange(len(m))]
		err2 = [err[i]/(m[1]*i) for i in xrange(len(m))]
		plt.errorbar(range(len(m)), m2, err2)
	else:
		m = msid(pos_segments)
		m2 = [m[i]/(m[1]*i) for i in xrange(len(m))]
		plt.plot(range(len(m)), m2)
	plt.show()
	

def write_msid(pos_segments, nom = 'msid.dat', error=True):
	if error: m, err = msid_and_err(pos_segments)
	else: m = msid(pos_segments)
	print "Debut de l'ecriture du fichier %s..." %(nom)
	if not nom in os.listdir(os.getcwd()): open(nom, 'w').close()
	with open(nom,'w') as f:
		if error:
			f.write("# n msid msid(n)/n*msid(1) err errprc\n\n")
			f.write("\n".join(["%i %f %f %f %f" %(i, m[i], m[i]/(m[1]*i), err[i], err[i]/(m[i]*i)) for i in xrange(1, len(m))]))
		else:
			f.write("# n msid msid(n)/n*msid(1)\n\n")
			f.write("\n".join([str(i)+' '+str(m[i])+' '+str(m[i]/(m[1]*i)) for i in xrange(1, len(m))]))
			

def write_seq_msid(seq, graph0, nom = 'msid.dat'):
	"""!!!afin d'accelerer le processus, les ruptures ne sont pas incluses!!!
	"""
	at_seg = get_segments(graph0, data=True)
	stime = [str(g.graph["timestep"]) for g in seq]
	print "Prise en compte de la periodicite de la boite..."
	pos_seg_seq = []
	for i in xrange(len(seq)):
		box = np.array([seq[i].graph['xhi']-seq[i].graph['xlo'], seq[i].graph['yhi']-seq[i].graph['ylo'], seq[i].graph['zhi']-seq[i].graph['zlo']])
		pos_seg_seq.append([])
		for j in xrange(len(at_seg)):
			pos_seg_seq[i].append([seq[i].node[at_seg[j][0]]['pos']])
			for k in xrange(1,len(at_seg[j])):
				pos_seg_seq[i][j].append(compenser_periode_box(pos_seg_seq[i][j][-1], seq[i].node[at_seg[j][k]]['pos'], box))
	print "Calcul du msid..."
	msid_seq = [msid(pos_seg) for pos_seg in pos_seg_seq]
	print "Debut de l'ecriture du fichier %s..." %(nom)
	if not nom in os.listdir(os.getcwd()): open(nom, 'w').close()
	with open(nom,'w') as f:
		f.write("#msid(1) pour tous t")
		f.write("\n1 "+" ".join([str(m[1]) for m in msid_seq]))
		f.write("\n#msid(n)/n*msid(1) pour tous t")
		f.write("\nn t"+' t'.join(stime))
		for n in xrange(1,max([len(s) for s in at_seg])):
			f.write("\n"+str(n)+" "+" ".join([str(m[n]/(m[1]*n)) for m in msid_seq]))
		
def vec_internes(pos_segments):
	for pos_section in pos_segments:
		for i in xrange(1,len(pos_section)):
			vec = pos_section[i] - pos_section[i-1]
			norm = np.linalg.norm(vec)
			if norm > 1.3: print norm, vec
		print 'bout a bout', np.linalg.norm(pos_section[-1]-pos_section[0])
		
def densities(dclt=None, dclm=None, dclh=None, Vt=None, Vm=None, Vh=None, rVmh=None, rdmh=None, f=None, c=None):
	"""densities(dclt=None, dclm=None, dclh=None, Vt=None, Vm=None, Vh=None, rVmh=None, rdmh=None, f=None, c=None)-> None
		(float, float, float, float, float, float, float, float, float, float)-> None
	Affiche les valeurs de dclt, dclm, dclh, Vt, Vm, Vh, rVmh, rdmh, f et c
	Doit etre rempli: Vt, et 3 parametres supplementaires dont un parametre supplementaire en volume et deux parametres en densite de reticulation  
	t:boite entiere, m:matrice, h:elements inclus dans la matrice
	Vx: volume de la boite occupe par x, dclx: densite de reticulation dans x, rVmh=Vm/Vh, rdmh=dclm/dclh, f = Vh/Vt, c = dclh/dclm
	"""
	if len([0 for i in (dclt, dclm, dclh, Vm, Vh, rVmh, rdmh, f, c) if i]) < 3 or not Vt: 
		print 'il manque une information'
		return
	
	#calcul des volumes a partir de l'info volume
	if f: Vh = f*Vt
	if rVmh: Vh = Vt/(1.+rVmh)
	if Vh: Vm = Vt-Vh
	elif Vm: Vh = Vt-Vm
	f = Vh/Vt
	rVmh = Vm/Vh
	#infos acquises Vt, Vh, Vm, f et rVms
	
	#calcul des densites a partir des 2 infos densites
	if rdmh: c = 1./rdmh
	if c:
		if dclt: dclm = dclt*Vt/(Vm+Vh*c)
		if dclm: dclh = c*dclm
		elif dclh: dclm = dclh/c
		dclt = (dclm*Vm + dclh*Vh)/Vt 
	elif dclt and dclm: dclh = (dclt*Vt - dclm*Vm)/Vh
	elif dclt and dclh: dclm = (dclt*Vt - dclh*Vh)/Vm
	elif dclm and dclh: dclt = (dclm*Vm + dclh*Vh)/Vt
	c = dclh/dclm
	rdmh = 1./c
	#infos acquises dclt, dclh, dclm, c et rdmh	
		
	print 'densite de reticulation totale:',dclt, '\nde la matrice:',dclm, "\nde l'heterogeneite:",dclh, '\nVolume total:',Vt, '\nde la matrice:',Vm, "\nde l'heterogeneite:",Vh, "\nfraction volumique occupee par l'heterogeneite:",Vh/Vt, '\ncontraste de reticulation:', 1./rdmh
	return {'dclt':dclt, 'dclm':dclm, 'dclh':dclh, 'Vt':Vt, 'Vm':Vm, 'Vh':Vh, 'rVmh':rVmh, 'rdmh': rdmh, 'f':Vh/Vt, 'c':1./rdmh}

if __name__ == '__main__':	
	#routine pour verifier si toutes les fonctions repondent
	l_moy_segments(debug_pos_secs)
	barycentres(debug_pos_secs)
	sq_dist_extrem_seg(debug_pos_secs)
	sq_rayon_giration(debug_pos_secs)
	msid(debug_pos_secs)
	msid(debug_pos_secs)
	write_msid(debug_pos_secs)
	vec_internes(debug_pos_secs)
	

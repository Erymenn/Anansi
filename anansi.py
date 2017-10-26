#!/usr/bin/python -i

"""
---------------------------------------------------------
                    Anansi 1.0 17/07
------------------------------------------------------------
Code for creation and caracterisation of elastomer networks
Coded by Morgane Mahaud for Mateis lab
Contact: morganemahaud@hotmail.com

Versions: python 2.7 networkx 1.11

Files: 
anansi.py
anansi_active_graph_fonctions.py
anansi_crosslink_particularities.py
anansi_file_conversion.py
anansi_math_fonctions.py
anansi_passive_graph_fonctions.py
debug_data
debug_dumpbd
------------------------------------------------------------
"""

#LIBRARIES AND FONCTIONS IMPORT
import time
import matplotlib.pyplot as plt
from math import sqrt
from anansi_active_graph_fonctions import *
from anansi_file_conversion import *
from anansi_math_fonctions import *
from anansi_passive_graph_fonctions import *


#DEFINITION OF FONCTIONS USED IN SCIPTS
def create_random_heterogeneous_crosslinking(graph, data_output, b, c, shape, dclm, k, f):#!!!
    """(nx.Graph, str, b, c, shape, dclm, k, f) -> None
    Create an heterogeneously crosslinked elastomer from a polymer 
    shapes of heterogeneites: cross or sphere"""
    d = densities(Vt=b[2][0]*b[2][1]*b[2][2], dclm=dclm, c=k, f=f)
    if shape == 'croix': z = croix(graph, c, b, f)
    elif shape == 'sphere': z = sphere(graph, c, b, f=f)
    else: print "heu", shape
    m = list(set(graph.nodes()).difference(set(z)))
    modified_graph = reticuler(graph, ['miniboucle'], d['dclh'], zone = z, visuel = False, blabla=True)
    twice_modified_graph = reticuler(modified_graph, ['alea'], d['dclm'], zone = m, visuel = False, sauvdata = data_output)


#EXEMPLES OF SCRIPTS
#-----------------------
"""TO USE ONE SCRIPT: 

IN A NEW FILE IN THE SAME FOLDER (OR DIFFERENT BUT BE CAREFUL OF PATHS), WROTE:
#!/usr/bin/python -i
import anansi
*SCRIPT*
THEN EXECUTE THE FILE
    
OR IN THE PRESENT, UNCOMMENT ONLY THE NEEDED SCRIPT AND EXECUTE ANANSI.PY
"""
    
if __name__ == '__main__':
    
    pass
    #CREATE AN HOMOGENEOUS ELASTOMER FROM A POLYMER DATA LAMMPS FILE
    #---------------------------
    #file ='/home/morgane/Documents/execution/restart to data/data-refroidi-1-'+ch
    #g = data_to_graph(file)
    #h = reticuler(g, ['alea'], dr, visuel = False, blabla = False, sauvdata = "data-1-dr"+str(dr)+"-"+ch)


    #CREATE SEVERAL HOMOGENEOUS ELASTOMERS FROM POLYMER DATA LAMMPS FILES
    #---------------------------
    #for ch in ["100-1000"]:#, "500-200", "2000-50", "10000-10"]:
        #for dr in [0.05, 0.1]:#[0.002*i for i in range(5, 21)]:
            #fich ='/home/morgane/Documents/execution/restart to data/data-refroidi-1-'+ch
            #g = data_to_graph(fich)
            #print 100000./volume(g)
            #h = reticuler(g, ['alea'], dr, visuel = False, blabla = False, sauvdata = "data-1-dr"+str(dr)+"-"+ch)
            ##print ch, dr, "fs", "nuc"
            ##print str(dr)+"\t"+str(fraction_soluble(h))+"\t"+str(get_densite_actives_chimiques(h))



    #CREATE AN ELASTOMER WITH SPHERIC HETEROGENEITY FROM A POLYMER DATA LAMMPS FILE

    #CREATE AN ELASTOMER WITH CROSS HETEROGENEITY FROM A POLYMER DATA LAMMPS FILE

    #CREATE A DOUBLE NETWORK FROM A POLYMER DATA LAMMPS FILE

    #CALCULATE 











    #fich = "/home/morgane/Documents/data/refroidi/data-refroidi-1-100-1000"
    #graph = data_to_graph(fich)
    #b = boxparam(graph)
    #c = cadrillage(graph, b)    
    #dcl = 0.02
    #for k in [2.5]:
        #for f in [0.2]:
            #for shape in ["sphere"]:#, "croix"]:
                #graph = data_to_graph(fich)
                #creer_reticulation_aleatoire_heterogene(graph,'data-'+shape+'2-dcl0.02-f'+str(f)+'k'+str(k), b, c, shape, dcl, k, f)
                #creer_reticulation_aleatoire_heterogene(graph,'test', b, c, shape, dcl, k, f)

    #calculer fraction soluble et decompter particularites_to_type
    #----------------------------
    #for k in [1, 2.5, 5]:
       #for f in [0.5]:
           #for shape in ["croix"]:#, "croix"]:
               #fich = 'data2-'+shape+'-dcl0.02-f'+str(f)+'k'+str(k)
               #print shape, f, k
               #g = data_to_graph(fich)
               #print get_densite_actives_chimiques(g)
               #h = get_graph_actif(g)
               #print get_densite_actives_chimiques(h)
               #decompter_particularites(h)
               
    #for name in ["data-db50-dr0.02-1000-100", "data-db50-dr0.02-1000-100-trac1", "data-db50-dra0.029-drb0.011-1000-100", "data-db50-dra0.029-drb0.011-1000-100-trac1", "data-db50-dra0.033-drb0.007-1000-100", "data-db50-dra0.033-drb0.007-1000-100-trac1"]:
       ##drc=0.04-dr
       #fich = "/home/morgane/Documents/data/rupture/db/data/"+name
       #g = data_to_graph(fich)
       #print name
       #print get_densite_actives_chimiques(g)
       #h = get_graph_actif(g)
       #decompter_particularites(h)

    #fich = "/home/morgane/Documents/data/rupture/data/data-dr0.04-100-1000"
    #g = data_to_graph(fich)
    #decompter_particularites(g)
    
    #for drc in [0.03]:
        #fich = "/home/morgane/Documents/data/refroidi/data-refroidi-1-100-1000"
        #nomsauv = "/home/morgane/Documents/data/rupture/data/data-dr0.04-drc"+str(drc)+"-100-1000.dat"
        ##fich = "/home/morgane/Documents/data/rupture/data/data-dr"+str(dr)+"-"+ch
        #graph = data_to_graph(fich)
        #h = reticuler(graph, ['alea'], 0.04, visuel = False, blabla = False)
        #g = couper(h, drc, sauvdata = nomsauv)
    
    #calculer distance end-end
    #-----------------------
    #for a in ["", "-trac1"]:
        #for b in ["dr0.02-1000-100", "dra0.029-drb0.011-1000-100", "dra0.033-drb0.007-1000-100"]:
            ##M = 100000/N
            #N=1000
            #drc=0
            #dr=0.02
            #fich="/home/morgane/Documents/data/rupture/db/data/data-db50-"+b+a
            ##fich = "/home/morgane/Documents/programmation/c_solver_polymer_md_mechanics/lib/graph/data-sphere-dclm0.02-f"+str(fv)+"k"+str(k)
            #g = data_to_graph(fich)
            #pos_sec = get_sections(g, pos = True)
            #ee = dist_extrem_sec(pos_sec)
            #Nsc = float(N)/(1+2*dr*N+drc*N)
            #print a, b, ee, Nsc, Nsc/ee
            #
            #nuc = get_densite_actives_chimiques(g)
            #h = get_graph_actif(g)
            #pos_sec = get_sections(h, pos = True)
            #eea = dist_extrem_sec(pos_sec)
            #print "tot", N, dr-drc, nuc, ee, sqrt(Nsc), Nsc/ee#, sqrt(Nsca), Nsca/ee
            
    #nom = "topo-sphere_dclm-varie.dat"
    #fich = "/home/morgane/Documents/data/refroidi/data-refroidi-1-100-1000"
    #graph = data_to_graph(fich)
    #b = boxparam(graph)
    #c = cadrillage(graph, b)
    #if not nom in os.listdir(os.getcwd()): open(nom, 'w').close()    
    #with open(nom,'w') as f:
        #f.write("dclm f k fs fp nuc Nsct Nsca Nscp Nscs\n")
        #for frac in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:    
            #for k in [1, 2.5, 5]:
                #dclm=0.02
                #h = data_to_graph(fich)
                ##fich="data-sphere-dclm0.02-f%sk%s" %(str(frac), str(k))
                ##fich = "/home/morgane/Documents/data/rupture/data/data-dr"+str(dr)+"-"+ch
                ##h=data_to_graph(fich)
                #creer_reticulation_aleatoire_heterogene(h, 'test', b, c, "sphere", dclm, k, frac)
                ##h = couper(graph, drc)
                ##h = reticuler(graph, ['alea'], dr, visuel = False, blabla = False)#, sauvdata = "/home/morgane/Documents/data/topo/data-dr"+str(dr)+"-"+ch)
                ##nom_topo = "/home/morgane/Documents/programmation/c_solver_polymer_md_mechanics/scripts/topo/data_topologie-"+ch+"-dr"+str(dr)+"-drc"+str(drc)+".dat"
                #nom_topo = "/home/morgane/Documents/programmation/c_solver_polymer_md_mechanics/scripts/topo/data_topologie-"+"sphere"+"-dclm"+str(dclm)+"-f"+str(f)+"-k"+str(k)+".dat"
                ##nom_topo ="topo-"+fich
                #fs, fp = marquer_pendants_solubles(h)
                #nuc = get_densite_actives_chimiques(h)
                #rep = repartition_sections_topologie(h, nom_topo)
                ##N = int(ch.split("-")[1])
                #Nsct = float(sum([i*rep[0][i] for i in range(len(rep[0]))]))/max(1, sum([i for i in rep[0]]))
                #Nsca = float(sum([i*rep[1][i] for i in range(len(rep[1]))]))/max(1, sum([i for i in rep[1]]))
                #Nscp = float(sum([i*rep[2][i] for i in range(len(rep[2]))]))/max(1, sum([i for i in rep[2]]))
                #Nscs = float(sum([i*rep[3][i] for i in range(len(rep[3]))]))/max(1, sum([i for i in rep[3]]))
                #print dclm
                #print "fraction soluble:", fs
                #print "fraction de chaines pendantes:", fp
                #print "densite de chaines actives:", nuc
                #print "longueur moyenne des chaines:"
                #print "\t-totales: ", Nsct
                #print "\t-actives: ", Nsca
                #print "\t-pendantes: ", Nscp
                #print "\t-solubles: ", Nscs
                #f.write(" ".join([str(i) for i in [dclm,frac, k, fs, fp, nuc, Nsct, Nsca, Nscp, Nscs]])+"\n")
            
    
    #nom = "topo-croix.dat"
    #if not nom in os.listdir(os.getcwd()): open(nom, 'w').close()
    #with open(nom,'w') as f:
        #f.write("N dr drc fs fp nuc Nsct Nsca Nscp Nscs\n")
        #for ch in ['100-1000']:    
            #for name in ["data-croix-dcl0.02-f0.2k2.5","data-croix-dcl0.02-f0.2k5","data-croix-dcl0.02-f0.5k2.5","data-croix-dcl0.02-f0.5k5"]:
                #dr = 0
                #drc=0.02
                #fich=name
                ##fich = "/home/morgane/Documents/data/refroidi/data-refroidi-1-100-1000"
                ###fich = "/home/morgane/Documents/data/rupture/data/data-dr"+str(dr)+"-"+ch
                ##graph = data_to_graph(fich)
                #h=data_to_graph(fich)
                ##h = couper(graph, drc)
                ##h = reticuler(graph, ['alea'], dr, visuel = False, blabla = False)#, sauvdata = "/home/morgane/Documents/data/topo/data-dr"+str(dr)+"-"+ch)
                ##nom_topo = "/home/morgane/Documents/programmation/c_solver_polymer_md_mechanics/scripts/topo/data_topologie-"+ch+"-dr"+str(dr)+"-drc"+str(drc)+".dat"
                #nom_topo ="topo-"+name
                #fs, fp = marquer_pendants_solubles(h)
                #nuc = get_densite_actives_chimiques(h)
                #rep = repartition_sections_topologie(h, nom_topo)
                ##N = int(ch.split("-")[1])
                #Nsct = float(sum([i*rep[0][i] for i in range(len(rep[0]))]))/max(1, sum([i for i in rep[0]]))
                #Nsca = float(sum([i*rep[1][i] for i in range(len(rep[1]))]))/max(1, sum([i for i in rep[1]]))
                #Nscp = float(sum([i*rep[2][i] for i in range(len(rep[2]))]))/max(1, sum([i for i in rep[2]]))
                #Nscs = float(sum([i*rep[3][i] for i in range(len(rep[3]))]))/max(1, sum([i for i in rep[3]]))
                #print ch, dr, drc
                #print "fraction soluble:", fs
                #print "fraction de chaines pendantes:", fp
                #print "densite de chaines actives:", nuc
                #print "longueur moyenne des chaines:"
                #print "\t-totales: ", Nsct
                #print "\t-actives: ", Nsca
                #print "\t-pendantes: ", Nscp
                #print "\t-solubles: ", Nscs
                #f.write(" ".join([str(i) for i in [ch, dr, drc, fs, fp, nuc, Nsct, Nsca, Nscp, Nscs]])+"\n")
                    
    #nom = "topo-sphere-influef-c1-N1000-dclt0.02.dat"
    #nom = "temp"
    #if not nom in os.listdir(os.getcwd()): open(nom, 'w').close()
    #with open(nom,'w') as fil:
        #fil.write("N drt f c drm drh fs fp nuc Nsct Nsca Nscp Nscs\n")
        ##for ch in ['100-1000']:
            ##for dr in [0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04]:
        #for f in [0.6]:
            #for c in [5]:
                #N = 200
                #dclt =0.02
                #d = densities(Vt=110600, dclt = dclt, f = f, c = c)
                ##fich = "/home/morgane/Documents/data/rupture/heterogene/data/data-sphere1-rVms"+str(ch[0])+"-drm"+str(ch[1])+"-drs"+str(ch[2])+"-T1-100-1000"
                #fich = "/home/morgane/Documents/data/refroidi/data-refroidi-1-500-200"
                #g = data_to_graph(fich)
                #b = boxparam(g)
                #cad = cadrillage(g, b)
                #z = croix(g, cad, b, d['f'])
                ##z = sphere(g, cad, b, d['f'])
                #z2 = list(set(g.nodes()).difference(set(z)))
                #h = reticuler(g, ['alea'], d['dclh'], zone = z, visuel = False, blabla = False)
                #graph = reticuler(h, ['alea'], d['dclm'], zone = z2, visuel = False, blabla = False)
                #nom_topo = "/home/morgane/Documents/programmation/c_solver_polymer_md_mechanics/scripts/topo/data_topologie-croix-f"+str(f)+"-c"+str(c)+"N"+str(N)+".dat"
                ##nom_topo = "/home/morgane/Documents/programmation/c_solver_polymer_md_mechanics/scripts/topo/test2.dat"
                #fs, fp = marquer_pendants_solubles(graph)
                #nuc = get_densite_actives_chimiques(graph)
                #rep = repartition_sections_topologie(graph, nom_topo)
                #Nsct = float(sum([i*rep[0][i] for i in range(len(rep[0]))]))/max(1, sum([i for i in rep[0]]))
                #Nsca = float(sum([i*rep[1][i] for i in range(len(rep[1]))]))/max(1, sum([i for i in rep[1]]))
                #Nscp = float(sum([i*rep[2][i] for i in range(len(rep[2]))]))/max(1, sum([i for i in rep[2]]))
                #Nscs = float(sum([i*rep[3][i] for i in range(len(rep[3]))]))/max(1, sum([i for i in rep[3]]))
                #print N, f, c
                #print "fraction soluble:", fs
                #print "fraction de chaines pendantes:", fp
                #print "densite de chaines actives:", nuc
                #print "longueur moyenne des chaines:"
                #print "\t-totales: ", Nsct
                #print "\t-actives: ", Nsca
                #print "\t-pendantes: ", Nscp
                #print "\t-solubles: ", Nscs
                #fil.write(" ".join([str(i) for i in [N, d['dclt'], d['f'], d['c'], d['dclm'], d['dclh'], fs, fp, nuc, Nsct, Nsca, Nscp, Nscs]])+"\n")
    
    #fich = '/home/morgane/Documents/data/rupture/heterogene/data/data-sphere5-alea-dr0.0175-0.0005-T1-500-200'
    #dumpbd = None
    #dumpid = 'dump-0.1-0.1.dat'
    #g = dump_to_seq(dumpid)[0]
    #print '\n'.join(compter_par_caracnode(g, ['type','mol']).values())
    
    #creer data 1 mol
    #---------------------
    #fich='/home/morgane/Documents/data/refroidi/data-refroidi-1-100-1000'
    #g=data_to_graph(fich)
    #m=1
    #for n in g.nodes():
        #if g.node[n]['mol'] != m: g.remove_node(n)
    #unwrap(g)
    #g.graph['atoms'] = g.order()
    #g.graph['bonds'] = g.size()
    #print g.graph
    #graph_to_data(g, "/home/morgane/Documents/data/refroidi/data-refroidi-1-100-1000_mol1")
    
    
    ###creer croix ou sphere 29/01/16
    ##--------------------
    #fich = '/home/morgane/Documents/data/refroidi/data-refroidi-1-100-1000'
    ##fich ='/home/morgane/Documents/programmation/reseaux/data-sphere1-rVms3.0-drm0.010-drs0.050-T1-100-1000'
    #nomsauv = '5spheres'
    ##nomsauv = "sphere1-rVms%.1f-drm%.3f-drs%.3f-T1-100-1000" %(rVms, drm, drs)
    #g = data_to_graph(fich)
    #b = boxparam(g)
    #c = cadrillage(g, b)
    ##d = densities(Vt=110590, drm=0.010, rVms = 3., drs = 0.05)
    ##drt = d['drt']
    ##drs = d['drs']
    ##drm = d['drm']
    ##rVms = d['rVms']
    ##rmax = 1.5
    ##-----------croix
    ##z = croix(g, c, b, rVms)
    ##-----------sphere
    ##z = sphere(g, c, b, rVms)
    ##-----------5spheres
    #z = sphere(g, c, b, centre = [7, 7, 7], rayon = 12.)
    #z.extend(sphere(g, c, b, centre = [21, 33, 10], rayon = 12.))
    #z.extend(sphere(g, c, b, centre = [20, 10, 30], rayon = 12.))
    #z.extend(sphere(g, c, b, centre = [40, 21, 21], rayon = 12.))
    #z.extend(sphere(g, c, b, centre = [4, 35, 33], rayon = 12.))
    ##------------------
    #s2 = list(set(g.nodes()).difference(set(z)))
    ###h = reticuler(g, ['alea'], drs, rmax, 3, zone = z, visuel = False)
    ###h = reticuler(g, ['alea'], drm, rmax, 3, zone = s2, visuel = False, sauvdata = 'data-'+nomsauv)
    #graph_to_dump(g, 'dump-'+nomsauv)
    #suivre_zones_from_dump('dump-'+nomsauv, (z, s2))
    
    #creer fissure
    #-------------------
    #fich = "/home/morgane/Documents/data/rupture/data/data-dr0.02-seuil3-1-100-1000"
    #nomsauv = "/home/morgane/Documents/data/rupture/data/data-fissure-dr0.02-seuil3-1-100-1000"
    #g = data_to_graph(fich)
    #b = boxparam(g)
    #c = cadrillage(g, b)
    #rVst = 0.005
    #boxhi, boxlo, box, l, inf = b
    #centre = np.array([boxlo[0] + box[0]/2, boxlo[1] + box[0]/2, boxlo[2]])
    #z = sphere(g, c, b, rVst = rVst, centre = centre)
    #retirer_zone(g, z)
    #graph_to_data(g, nomsauv)
    #graph_to_dump(g, 'dump-fiss')
    
    
    #suivre rupture heterogene
    #-------------------
    #dump = '/home/morgane/Documents/data/rupture/heterogene/croix/croix-rVms3.0-drm0.010-drs0.050-T1-100-1000/dump-croix-rVms3.0-drm0.010-drs0.050-T1-100-1000bis.dat'
    #dumpbd = '/home/morgane/Documents/data/rupture/heterogene/croix/croix-rVms3.0-drm0.010-drs0.050-T1-100-1000/dump-bd-croix-rVms3.0-drm0.010-drs0.050-T1-100-1000bis.dat'
    #data = '/home/morgane/Documents/data/rupture/heterogene/data/data-croix-rVms3.0-drm0.010-drs0.050-T1-100-1000'
    #data = '/home/morgane/Documents/data/rupture/heterogene/data/data-sphere1-rVms3.0-drm0.010-drs0.050-T1-100-1000'
    #dump = ['/home/morgane/Documents/data/rupture/heterogene/sphere/sphere1-rVms3.0-drm0.010-drs0.050-T1-100-1000/dump-sphere1-rVms3.0-drm0.010-drs0.050-T1-100-1000.dat',\
        #'/home/morgane/Documents/data/rupture/heterogene/sphere/sphere1-rVms3.0-drm0.010-drs0.050-T1-100-1000/dump-sphere1-rVms3.0-drm0.010-drs0.050-T1-100-1000bis.dat',\
        #'/home/morgane/Documents/data/rupture/heterogene/sphere/sphere1-rVms3.0-drm0.010-drs0.050-T1-100-1000/dump-sphere1-rVms3.0-drm0.010-drs0.050-T1-100-1000ter.dat']
    #dumpbd = ['/home/morgane/Documents/data/rupture/heterogene/sphere/sphere1-rVms3.0-drm0.010-drs0.050-T1-100-1000/dump-bd-sphere1-rVms3.0-drm0.010-drs0.050-T1-100-1000.dat',\
           #'/home/morgane/Documents/data/rupture/heterogene/sphere/sphere1-rVms3.0-drm0.010-drs0.050-T1-100-1000/dump-bd-sphere1-rVms3.0-drm0.010-drs0.050-T1-100-1000bis.dat',\
           #'/home/morgane/Documents/data/rupture/heterogene/sphere/sphere1-rVms3.0-drm0.010-drs0.050-T1-100-1000/dump-bd-sphere1-rVms3.0-drm0.010-drs0.050-T1-100-1000ter.dat']
    #nomsauv = "rupture_sphere1-rVms3.0-drm0.010-drs0.050-T1-100-1000"
    #g = data_to_graph(data)
    #b = boxparam(g)
    #c = cadrillage(g, b)
    #rVms = 3
    #rmax = 1.5
    ##z = croix(g, c, b, rVms)
    #z = sphere(g, c, b, rVms)
    #s2 = list(set(g.nodes()).difference(set(z)))
    #suivre_rupture_from_dumps(dump, dumpbd, nomsauv)
    #suivre_zones_from_dump(nomsauv, (z, s2), nomsauv+"-zones")
    
    #suivre rupture
    #-----------------
    #dumpbd = ["/home/morgane/Documents/data/rupture/rF100_hetero/sphere2-dcl0.02-f0.2k5/dump-bd-sphere2-dcl0.02-f0.2k5.dat"]
    #dump = ["/home/morgane/Documents/data/rupture/rF100_hetero/sphere2-dcl0.02-f0.2k5/dump-sphere2-dcl0.02-f0.2k5.dat"]
    #nomsauv = "dumprup-sphere2-dcl0.02-f0.2k5"
    #fich = "/home/morgane/Documents/data/refroidi/data-refroidi-1-100-1000"
    #g = data_to_graph(fich)
    #b = boxparam(g)
    #c = cadrillage(g, b)
    #d = densities(Vt=b[2][0]*b[2][1]*b[2][2], dclt=0.02, c=5, f=0.2)
    #z = sphere(g, c, b, 0.2)
    ##z = croix(g, c, b, 0.5)
    #m = list(set(g.nodes()).difference(set(z)))
    #suivre_rupture_from_dumps(dump, dumpbd, nomsauv)
    #suivre_zones_from_dump(nomsauv, (z, m), nomsauv+"-zones")
    
    #suivre rupture reseau db
    #-----------------
    #dumpbd = ["/home/morgane/Documents/data/rupture/db/db50-dra0.033-drb0-1000-100/dump-bd-db50-dra0.033-drb0-1000-100.dat","/home/morgane/Documents/data/rupture/db/db50-dra0.033-drb0.007-1000-100-trac1/dump-bd-db50-dra0.033-drb0.007-1000-100-trac1.dat", "/home/morgane/Documents/data/rupture/db/db50-dra0.033-drb0.007-1000-100-trac1/dump-bd-db50-dra0.033-drb0.007-1000-100-trac1bis.dat"]
    #dump = ["/home/morgane/Documents/data/rupture/db/db50-dra0.033-drb0-1000-100/dump-db50-dra0.033-drb0-1000-100.dat","/home/morgane/Documents/data/rupture/db/db50-dra0.033-drb0.007-1000-100-trac1/dump-db50-dra0.033-drb0.007-1000-100-trac1.dat", "/home/morgane/Documents/data/rupture/db/db50-dra0.033-drb0.007-1000-100-trac1/dump-db50-dra0.033-drb0.007-1000-100-trac1bis.dat"]
    #nomsauv = "dumprup-db50-dcl0.02-k5-trac"
    #suivre_rupture_from_dumps(dump, dumpbd, nomsauv)
    
    #mesure echevetrement
    #-----------------
    #data = "/home/morgane/Documents/Polymer-500-200-timestep465593.dat"
    ##dump = "/home/morgane/P2chpd/mmahaud/execution/ppa/data-refroidi-0.1-500-200/dump.dat"
    ##g = dump_to_seq(dump)[23]
    #g = data_to_graph(data)
    #pos_sections = get_sections(g, pos = True)
    #echant_lpp=3
    #lpp = lpp(pos_sections, echant_lpp)
    #moyn = moy_n(pos_sections)
    #sqendend = sq_end_end(pos_sections)
    #print "lpp:", lpp, "moyn:", moyn, "sqendend:", sqendend, "Ne:", moyn*sqendend/lpp**2#moyn*200/(6*lpp**2)#
    
    #reseau double
    #--------------------
    #data = "/home/morgane/Documents/execution/restart to data/data-db50-dra0.033-drb0-1000-100-trac1"
    #g= data_to_graph(data)
    #z1 = []
    #z2 = []
    #for n in g:
        #if g.node[n]['mol'] < 51: z1.append(n)
        #else: z2.append(n)
    #dr = 0.02 
    #c = 5
    #dra = 2*c*dr/(1+c)
    #drb = 2*dr/(1+c)
    #nomsauv="db50-dra%.3f-drb%.3f-1000-100" %(dra, drb)
    #h = reticuler(g, ['alea'], dra, 1.5, 3, zone = z1, visuel = False, sauvdata = 'data-'+nomsauv)
    #h = reticuler(h, ['alea'], drb, 1.5, 3, zone = z2, visuel = False, bondtype = 3, sauvdata = 'data-'+nomsauv)
    
    
    #suivre pendants solubles
    #---------------------------
    #for ch in ["10000-10", "500-200","100-1000"]:
    #for ch in ["100-1000"]:
        #for dr in [0.005, 0.01, 0.04]:
            #fich = "/home/morgane/Documents/data/rupture/data/data-dr"+str(dr)+"-"+ch
            #fich="/home/morgane/Documents/data/rupture/data/data-dr0.02-seuil3-1-"+ch
            #dump = "/home/morgane/Documents/data/rupture/dumprF100/dump-dr"+str(dr)+"-"+ch+".dat"
            #g = data_to_graph(fich)
            #marquer_pendants_solubles(g, dump)
            #print "fraction sol:", fraction_soluble(g)
            
    #longueur moyenne des sections de graph
    #---------------------------
    #for ch in ["500-200","100-1000"]:
        #for dr in [0.005, 0.01, 0.02, 0.04]:
            #fich = "/home/morgane/Documents/data/rupture/data/data-dr"+str(dr)+"-"+ch
            #g = data_to_graph(fich)
            ##print get_sections(get_graph_actif(g))
            #N = float(ch.split("-")[1])
            #M = float(ch.split("-")[0])
            #Nmoy = N/(dr*200000/M)
            #print ch, dr, "Nmoy_tot:", N_moy_sections(g), "Nmoy:", Nmoy #"Nmoy_act:", N_moy_sections(get_graph_actif(g)) 
    
                
            ##decompter_particularites(h)
            ##print "soluble:", fraction_soluble(g)
            ##print "pendante:", fraction_pendante(g)
    
    #calculer fraction soluble et decompter particularites_to_type
    #----------------------------
    #nom = "topo-db.dat"
    #if not nom in os.listdir(os.getcwd()): open(nom, 'w').close()
    #with open(nom,'w') as f:
        #f.write("fs fp nuc Nsct Nsca Nscp Nscs\n")
        #for t in ["", "-trac1"]:
            ###for dr in [0.0005*i for i in range(1,10)]:
            #for n in ["db50-dr0.02-1000-100", "db50-dra0.029-drb0.011-1000-100", "db50-dra0.033-drb0.007-1000-100"]:
                #fich = "/home/morgane/Documents/data/rupture/db/data/data-"+n+t
                #nom_topo = 'topo-'+n+t
                ##fich ="/home/morgane/Documents/data/topo/"+ch+"/data-dr"+str(dr)+"-"+ch
                #h = data_to_graph(fich)
                ##decompter_particularites(g)
                #print n+t
                ##fs, fp = marquer_pendants_solubles(h, db=True)
                #nuc = get_densite_actives_chimiques(h)
                #rep = repartition_sections_topologie(h, nom_topo)
                #Nsct = float(sum([i*rep[0][i] for i in range(len(rep[0]))]))/max(1, sum([i for i in rep[0]]))
                #Nsca = float(sum([i*rep[1][i] for i in range(len(rep[1]))]))/max(1, sum([i for i in rep[1]]))
                #Nscp = float(sum([i*rep[2][i] for i in range(len(rep[2]))]))/max(1, sum([i for i in rep[2]]))
                #Nscs = float(sum([i*rep[3][i] for i in range(len(rep[3]))]))/max(1, sum([i for i in rep[3]]))
                #print "fraction soluble:", fs
                #print "fraction de chaines pendantes:", fp
                #print "densite de chaines actives:", nuc
                #print "longueur moyenne des chaines:"
                #print "\t-totales: ", Nsct
                #print "\t-actives: ", Nsca
                #print "\t-pendantes: ", Nscp
                #print "\t-solubles: ", Nscs
                #f.write(" ".join([str(i) for i in [fs, fp, nuc, Nsct, Nsca, Nscp, Nscs]])+"\n")
                
    #for t in ["", "-trac1"]:
        #for n in ["db50-dr0.02-1000-100", "db50-dra0.029-drb0.011-1000-100", "db50-dra0.033-drb0.007-1000-100"]:
            #fich = "/home/morgane/Documents/data/rupture/db/data/data-"+n+t
            #h = data_to_graph(fich)
            #print n+t
            #nuc = get_densite_actives_chimiques(h, db=True)
            #print nuc
            
    #for dcl in [0.002*i for i in range(1,26)]:
        #fich = "/home/morgane/Documents/data/refroidi/data-refroidi-1-100-1000"
        #g = data_to_graph(fich)
        #h = reticuler(g, ['alea'], dcl, visuel = False, blabla = False)
        #print "dcl=",dcl
        #decompter_particularites(h, seuil=3)
            
    
    #fich = '/home/morgane/Documents/data/RLP/benchmark/M500N200-cutoff2.5/slurm-97907.out'
    #slurm_to_suivi_retic(fich, 'suivi_retic_R2.5', 100000) 
    
    #fich = "/home/morgane/Documents/execution/SOMM_LINUX/HotMono_new/Liquid_48M.dat"
    #somm_to_data(fich, 'data-RLP-bath48M')
    
    
    #dp = "/home/morgane/Documents/data/RLP/equilibre/rlp-10000-10/test.dat"
    #chain = "500-200"
    #dp = "/home/morgane/Documents/data/RLP/equilibre/rlp-"+chain+"/dump-rlp-"+chain+".dat"
    #tmp = dump_to_seq(dp)
    #seq = [tmp[i] for i in range(0,len(tmp),10)]
    #seq_to_dump(seq, "/home/morgane/Documents/data/RLP/equilibre/rlp-"+chain+"/dump-div10-rlp-"+chain+".dat")
    ##data0 = "/home/morgane/Documents/data/data-eq.dr0.02-seuil3-1-2000-50"
    
    #calculer un msid et son erreur
    ##-----------------
    #chain = "500-200"
    ##dp = "/home/morgane/Documents/data/RLP/equilibre/rlp-10000-10/dump-rlp-10000-10prem.dat"
    ##dp = "/home/morgane/Documents/data/RLP/equilibre/rlp-10000-10/dump-div10-rlp-10000-10.dat"
    ##dp = "/home/morgane/Documents/data/RLP/equilibre/rlp-2000-50/dump-t500000.dat"
    ##dp = "/home/morgane/Documents/data/RLP/N1000M200/restart.100000000"
    ##dp = "/home/morgane/Documents/data/RLP/equilibre/rlp-100-1000ter/dump-rlp-100-1000ter.dat"
    #data0 = "/home/morgane/Documents/execution/restart to data/data-200-1000"
    #g = data_to_graph(data0)
    ##seq = dump_to_seq(dp)[6]
    ##seq.add_edges_from(g.edges())
    #pos_sections = get_sections(g, pos=True)
    ##repartition_err(pos_sections)
    #write_msid(pos_sections, 'err'+chain+'.txt', True)
    ##see_msid(pos_sections, True)
    
    #calculer une serie de msid
    ##-----------------
    #chain = "100-1000"
    #dp = "/home/morgane/Documents/data/RLP/equilibre/rlp-100-1000qua/dump-rlp-100-1000qua.dat"
    #data0 = "/home/morgane/Documents/data/RLP/equilibre/rlp-"+chain+"/RLP-avant-equilibration-rlp-"+chain+".dat"
    #g = data_to_graph(data0)
    #seq = dump_to_seq(dp)
    #write_seq_msid(seq, g, "msid-rlp-"+chain+"qua.dat")
    
    #pos = get_sections(g, pos=True)
    #write_msid(pos, 'testmsid.txt')
    #for chain in ['10000-10']:
    ###for chain in ['500-200']:
        ##fich = '/home/morgane/Documents/execution/restart to data/data-refroidi-2-'+chain
        #fich = '/home/morgane/Documents/data/RLP/equilibre/eq-100-1000-Tp5000-T1.000/data-eq-Tp5000-T1.000-100-1000.80000000'
        #g = data_to_graph(fich)
        #pos = get_sections(g, pos=True)
        #write_msid(pos, 'msid80000000.dat')
        ##msid = msid2(pos)
        ##print 'debut msidh'
        ##x = range(1,len(msid))
        ##msidh = [msid[i]/(msid[1]*i) for i in x]
        ##plt.figure()
        ##plt.semilogx(x, msidh)
        ##plt.show()
        
    #calculer Rg
    #-----------------------
    #for ch in ['100-1000']:
        #fich = '/home/morgane/Documents/data/refroidi/data-refroidi-1-'+ch
        #g = data_to_graph(fich)
        #pos_sec = get_sections(g, pos = True)
        #sqrg = sq_rayon_giration(pos_sec)
        #sqd = dist_extrem_sec(pos_sec)**2
        #print sqrg, sqd, sqrg/sqd
    
    #fich = "/home/morgane/Documents/data/elastomere/trate/programme1/miniboucle0.03133/dump-miniboucle-0.5-prc0.03133-seuil3-500-200.dat"
    #datainit = "/home/morgane/Documents/data/elastomere/trate/programme1/miniboucle0.03133/data-miniboucle-0.5-prc0.03133-seuil3-500-200"
    #suivre_rupture_from_data(fich, datainit, 'test_suivi_rupture')
    
    #dps = ['/home/morgane/Documents/data/rupture/dump%s-alea-2-dr0.02-seuil3.dat' %(i) for i in ['', '-bd']]
    #t0 = time.clock()
    ##print dumpbd_to_bkbd_sets(dps[1])
    #suivre_particularites_from_dumps('test_suivi_rupture_dumps', dps[1], ['miniboucle'], 'test_suivi_part_dump')
    
    #suivi rupture et miniboucle
    #dump = '/home/morgane/Documents/data/rupture/alea-dr002-T1-100-1000-repet/dump-alea-dr002-T1-100-1000-repet-bis.dat'
    #dumpbd = '/home/morgane/Documents/data/rupture/alea-dr002-T1-100-1000-repet/dump-bd-alea-dr002-T1-100-1000-repet-bis.dat'
    ##dumpbd_to_weighted_list(dumpbd)
    #suivre_rupture_from_dumps(dump, dumpbd, nomsortie = 'rupture_dr0.02-1-100-1000-E590bis.dat')
    #suivre_particularites_from_dumps('sortie_rupture', dumpbd, ['miniboucle'], nomsortie = 'sortie_rupture_miniboucle', bd2eqclk = True, rompu = True)
    
    #calcul part fene/lj
    #-----------------------
    #dumpbd = '/home/morgane/Documents/data/rupture/dr0.02-1-100-1000/dump-bd-dr0.02-1-100-1000.dat'
    #seq = dumpbd_to_weighted_list(dumpbd, True)
    
    #fich = '/home/morgane/Documents/data/refroidi/data-refroidi-1-100-1000'
    #g = data_to_graph(fich)
    #dr = 0.02
    #rmax = 1.5
    #h = reticuler(g, ['alea'], dr, rmax, 3, visuel=False, sauvdata="F590-essai2-dr0.021-100-1000")
    #graph_to_dump(h, 'test')
    
    #1 sphere de rayon 20
    #fich = 
    #g = data_to_graph(fich)
    #b = boxparam(g)
    #c = cadrillage(g, b)
    #dr = 0.02
    #rmax = 1.5
    #s = sphere((b[1]+b[0])/2., 24, g, c, b)
    #h = reticuler(g, ['alea'], dr, rmax, 3, zone = s, visuel = False)
    #s2 = list(set(h.nodes()).difference(set(s)))
    #h = reticuler(g, ['alea'], dr, rmax, 3, zone = s2, visuel = False, sauvdata = 'sphere1-alea-dr0.01-0.01-T2-500-200')
    #graph_to_dump(h, 'test')
    
    #5 spheres de rayon 12 reparties aleatoirement
    #fich = '/home/morgane/Documents/execution/restart to data/data-refroidi-1-500-200'
    #g = data_to_graph(fich)
    #b = boxparam(g)
    #c = cadrillage(g, b)
    ##dr = 0.01
    ##drs = 0.002
    #dr = 0.01375
    #drs = 0.0004125
    #rmax = 1.5
    #s = sphere(7*np.ones(3), 12, g, c, b)
    #h = reticuler(g, ['alea'], drs, rmax, 3, zone = s, visuel = False)
    #s = sphere(np.array([21, 33, 10]), 12, g, c, b)
    #h = reticuler(h, ['alea'], drs, rmax, 3, zone = s, visuel = False)
    #s = sphere(np.array([20, 10, 30]), 12, g, c, b)
    #h = reticuler(h, ['alea'], drs, rmax, 3, zone = s, visuel = False)
    #s = sphere(np.array([40, 21, 21]), 12, g, c, b)
    #h = reticuler(h, ['alea'], drs, rmax, 3, zone = s, visuel = False)
    #s = sphere(np.array([4, 35, 33]), 12, g, c, b)
    #h = reticuler(h, ['alea'], drs, rmax, 3, zone = s, visuel = False)
    #s2 = list(set(h.nodes()).difference(set(s)))
    #h = reticuler(h, ['alea'], dr, rmax, 3, zone = s2, visuel = False, sauvdata = 'sphere5-alea-dr'+str(dr)+'-'+str(drs)+'-T1-500-200')
    #graph_to_dump(h, 'dump-sphere5-alea-dr'+str(dr)+'-'+str(drs)+'-T1-500-200')
    
    ##5 spheres de rayon 12 reparties aleatoirement
    #fich = '/home/morgane/Documents/execution/restart to data/data-refroidi-1-500-200'
    #dumpid = '/home/morgane/Documents/programmation/reseaux/rupture_sphere5-alea-dr0.0175-0.0005-T1-500-200'
    #g = data_to_graph(fich)
    #b = boxparam(g)
    #c = cadrillage(g, b)
    ##dr = 0.0175
    ##drs = 0.0005
    ##rmax = 1.5
    #z = []
    #z.append(sphere(7*np.ones(3), 12, g, c, b))
    #z.append(sphere(np.array([21, 33, 10]), 12, g, c, b))
    #z.append(sphere(np.array([20, 10, 30]), 12, g, c, b))
    #z.append(sphere(np.array([40, 21, 21]), 12, g, c, b))
    #z.append(sphere(np.array([4, 35, 33]), 12, g, c, b))
    #suivre_zones_from_dump(dumpid, z)
    
    ##fich = '/home/morgane/Documents/data/RLP/Polymer-500-200/Polymer-RLP-500-200-30000000.dat'
    #fich = 'data-refroidi-2-500-200'
    #fich = "data-alea-2-dr0.02-seuil3"
    ##fich = "/home/morgane/Documents/programmation/reseaux/Polymer-2000-50.dat"
    ##fich = "Polymer128000.dat"
    #g = data_to_graph(fich)
    ##print len(g)
    ##g.remove_nodes_from(range(100001, len(g)+1))
    ##print len(g)
    ##graph_to_data(g, "Polymer100000.dat")
    
    
    
    #particularites_to_type(g, ["miniboucle"], bonds = False, seuil = 3, bd2eqclk = True)
    #graph_to_dump(g, "dump-init-alea-2-dr0.02-seuil3")
    
    #print [n for n in g.nodes() if g.degree(n)>3]
    #s = get_sections(g, pos = True)
    #print s
    #print "msid1"
    #t1 = time.clock()
    #msid(s)
    #t2 = time.clock()
    #print "msid2"
    #msid2(s)
    #t3 = time.clock()
    #print "1:", t2-t1, "2:", t3-t2
    
    #boxhi = np.array([g.graph['xhi'], g.graph['yhi'], g.graph['zhi']])
    #boxlo = np.array([g.graph['xlo'], g.graph['ylo'], g.graph['zlo']])
    #box = boxhi-boxlo
    #sbis = get_sections(g, data = True)
    
    #for num in range(len(s)):
        #if g.node[sbis[num][0]]['mol'] == 10:
            #print "\n".join([str(i) for i in s[num]]), num
            #s3 = np.array([g.node[sbis[num][i]]['pos'] for i in range(len(sbis[num]))])
            #s1 = np.array(s[num])
        #if (s1-s3).any(): 
        #print s1-s3exit
    #vec_internes(s)
    #write_msid(s)
    #l = len (s)
    #print l_moy_segments(s)
    #print sum(barycentres(s))/l
    #print sq_dist_extrem_sec(s)
    #print sq_rayon_giration(s)
    #print 'debut msid'
    #msid = msid2(s)
    #write_msid(s)
    
    #for dr in [0.005, 0.01, 0.04]:
        #tag = '1-100-1000'
        #fich = '/home/morgane/Documents/data/refroidi/data-refroidi-'+tag
        #g = data_to_graph(fich)
        ##dr = 0.02
        #seuil = 3
        #reticuler(g, ['alea'], dr, 1.5, seuil = seuil, visuel = False, sauvdata = "./retic/data-dr"+str(dr)+"-seuil"+str(seuil)+'-'+tag)[1]
    #decompter_particularites(g)
    #particularites_to_type(g, ['intra', 'miniboucle', 'hexa'])
    #graph_to_data(g, "test")
    #dr = 0.04
    #reticuler(g, ['alea'], dr, 1.5, seuil = seuil, visuel = False, sauvdata = "data-*-2-dr"+str(dr)+"-seuil"+str(seuil))[1]

    #rep = repartition_sections(g)
    #print len(rep), rep[-1]
    #g = data_to_graph(fich)
    #nodelist = [i for i in g.nodes() if g.node[i]['mol']<10]
    #edgelist = []

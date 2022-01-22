from collections import defaultdict
from math import sqrt
import re

dictGenreNomLocus=defaultdict(list)
i=0
with open("Genre_Infect.tsv","r") as f3:
	for li in f3:
		if not li.startswith("Locus"):
			li=li.rstrip("\n")
			lp=li.split("\t")
			dictGenreNomLocus[lp[-4]].append(lp[0]) #Genre avec les noms de locus
			i+=1


dictProteinesNomVirus=defaultdict(str)
with open("virSorter_proteins_entete.faa", "r") as f2:
	for li in f2:
		li=li.rstrip("\n")
		lp=li.split()
		nomVirus=lp[1]
		proteine=lp[0].lstrip(">")
		dictProteinesNomVirus[proteine]=nomVirus #clef sont les protéines et en valeur l'id de virus associé à la protéine



dictNumeroGrpOrthoNomVirus=defaultdict(set)
i=0
with open("phrogs_list","r") as f1:
	for li in f1:
		li=li.rstrip("\n")
		i+=1
		lp=li.split()
		for virus in lp:
			match=re.search("_p\d+$",virus) #recherche des protéines finissant par _pEEEE
			if match:
				virus=virus.rstrip(match.group()) #si c'est une protéine finissant par _pEEEE, on enleve cette partie de l'identifiant et on ajoute au dictionnaire
				dictNumeroGrpOrthoNomVirus["phrogs_"+str(i)].add(virus)
			else:
				dictNumeroGrpOrthoNomVirus["phrogs_"+str(i)].add(dictProteinesNomVirus[virus]) #sinon on va chercher dans le dictionnaire créé avant le nom de virus associé à la protéine (VirSorter)


dictGenreNumGrpOrthoApparition=defaultdict(lambda:defaultdict(int))
for genre in dictGenreNomLocus:
	for nomLocus in dictGenreNomLocus[genre]:
		for numeroGrpOrtho in dictNumeroGrpOrthoNomVirus:
			if nomLocus in dictNumeroGrpOrthoNomVirus[numeroGrpOrtho]:
				dictGenreNumGrpOrthoApparition[genre][numeroGrpOrtho]+=1


i=0
for genre in dictGenreNumGrpOrthoApparition:
	for numeroGrpOrtho in dictGenreNumGrpOrthoApparition[genre]:
		 if i:
		 	with open("Conservation_grp_ortho.tsv","a") as fichier:
		 		fichier.write(genre + "\t" + numeroGrpOrtho + "\t" + str(dictGenreNumGrpOrthoApparition[genre][numeroGrpOrtho]/len(dictGenreNomLocus[genre])) + "\n")
		 else:
		 	i+=1
		 	with open("Conservation_grp_ortho.tsv","w") as fichier:
		 		fichier.write("Genre\tPhrogs\tConservation\n")


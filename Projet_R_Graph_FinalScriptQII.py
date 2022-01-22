from collections import defaultdict
from math import sqrt
import re



list_ortho = "phrogs_list"
prot_virus = "virSorter_proteins_entete.faa"
Grp_virus = "Grp_virus.tsv"

dict_idProt_ortho=defaultdict(str)
dict_idVirus_setOrtho=defaultdict(set)
dict_idProt_idVirus=defaultdict(str)


with open(list_ortho, "r") as f1 :
    compt=0
    for li in f1 :
        li=li.rstrip("\n")
        ls=li.split()
        compt+=1
        ortho="Ortho"+str(compt)
        ortho = compt   # numéro de ligne du groupe orthologue

        for prot in ls :
        	
            if prot.startswith("p") :
                dict_idProt_ortho[prot]=ortho
                # print(prot)
            else :
                id_virus = re.sub("_p\d+","",prot)
                #print(prot, id_virus)
                dict_idVirus_setOrtho[id_virus].add(ortho)


with open(prot_virus, "r") as f2 :
    for li in f2 :
        li = li.rstrip("\n")
        li = li.lstrip(">")
        ls = li.split()
        id_prot = ls[0]
        id_virus= ls[1]
        dict_idProt_idVirus[id_prot] = id_virus

for idProt in dict_idProt_ortho :
    ortho = dict_idProt_ortho[idProt]
    id_virus = dict_idProt_idVirus[idProt]
    dict_idVirus_setOrtho[id_virus].add(ortho)



dict_GroupeIdVirus_ortho = defaultdict(list)
dict_virus1_virus2_nbOrthoCommun = defaultdict(lambda : defaultdict(int))

with open("Grp_virus.tsv", "r") as f3 :
        for li in f3 :
            li=li.rstrip("\n")
            ls=li.split("\t")
            idVirus=ls[0]
            ortho = list(dict_idVirus_setOrtho[idVirus])
            dict_GroupeIdVirus_ortho[idVirus]=ortho

listDejaVue=[]

for virus1 in dict_GroupeIdVirus_ortho :
    listDejaVue.append(virus1)
    for virus2 in dict_GroupeIdVirus_ortho :
        if virus2 not in listDejaVue :
            for ortho1 in dict_GroupeIdVirus_ortho[virus1] :
                if ortho1 in dict_GroupeIdVirus_ortho[virus2] :
                    dict_virus1_virus2_nbOrthoCommun[virus1][virus2]+=1


repartition_grps_orthologues=open('repartition_grps_orthologues.tsv','w')
repartition_grps_orthologues.write("virus1"+"\t"+"virus2"+"\t"+ "nb_commun" +"\n" )
for virus1 in dict_virus1_virus2_nbOrthoCommun :
        for virus2 in dict_virus1_virus2_nbOrthoCommun[virus1] :
            nb_commun=dict_virus1_virus2_nbOrthoCommun[virus1][virus2]
            repartition_grps_orthologues.write(virus1+"\t"+virus2+"\t"+ str(nb_commun)+"\n" )

repartition_grps_orthologues.close()

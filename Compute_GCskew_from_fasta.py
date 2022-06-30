import os 
import matplotlib.pyplot as plt
import numpy as np

def GC_skew_SW(fasta_file, fenetre, step):
    """
    Parameters : 
    - fasta_file : (str) path et nom du fichier fasta à trouver
    - fenetre : (int) la taille de la fenêtre glissante permettant le calcul du GC_skew
    
    Returns : 
    - GC_list : (list) une liste contenant le GC_contenet pour chaque position dans la séquence
    - taille_seq : (int) qui est la taille de la séquence contenue dans le fichier fasta (optionnel)
    """
    pfile = open(os.getcwd()+fasta_file, 'r')
   
    #I - Initialisation des paramètres
    #GC_list = []
    C,G = 0,0
    content = ""
    t = 0
    #II - On extrait uniquement la séquence sans prendre en compte la première ligne 
    for line in pfile :
        if line[0]==">":
            head = line
        break
        
    content = pfile.read().replace(head,"").replace("\n","")
    
    #III - Calcul du GC_skew
    cumul = []
    taille_seq = len(content)
    p = 0
    GC_list = np.zeros(taille_seq)
    for nt in range(int(fenetre/2), taille_seq-int(fenetre/2),step):
        for i in range(fenetre):
            if content[nt-int(fenetre/2)+i]=="G":
                G += 1
            if content[nt-int(fenetre/2)+i]=="C":
                C +=1
        if (C+G)==0:
            GC_list[nt]=0
        else :
            GC_list[nt]=(C-G)/(C+G)    
        C,G = 0,0
        
    for i in range(len(GC_list)):
        if GC_list[i]!=0:
            t = GC_list[i]
        if GC_list[i]==0:
            GC_list[i]=t
            
    for i in range(len(GC_list)-1):
        cumul.append(GC_list[i]+GC_list[i+1])
            
    return GC_list,cumul

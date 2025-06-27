import matplotlib.pyplot as plt



def histogram(data):
    # Fusion des données pour avoir une seule entrée par clé
    merged = {}
    for d in data:
        for k, v in d.items():
            if k not in merged:
                merged[k] = []
            merged[k].extend(v)

    # Trier les clés
    sorted_keys = sorted(merged.keys())
    longueurs = [len(merged[k]) for k in sorted_keys]

    # Tracer l'histogramme (diagramme en barres)
    plt.figure(figsize=(10, 6))
    plt.bar(sorted_keys, longueurs, width=0.5, color='skyblue', edgecolor='black')
    plt.xlabel('Clés')
    plt.ylabel('Longueur des listes')
    plt.title('Histogramme des longueurs de listes par clé')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.show()
    
    
def histo(data):
    key = []
    val = []
    for dico in data:
        for k, v in dico.items():
            key.append(k)
            val.append(len(v))
    print(key)
    print(val)
    # Tracer l'histogramme (diagramme en barres)
    plt.hist(key, val)
    plt.plot()


    
if __name__ =='__main__':
    
    
    # Exemple de données
    data = [
        {1 : [1, 2, 3]},
        {2 : [4, 5]},
        {3 : [6]},
        {4.5 : [7, 8, 9, 10]},
        {5 : []},
        {6 : [11, 12]},
        {7.1 : [13, 14, 15, 16, 17]},
        {8.3 : [18]},
        {9.5 : [19, 20, 21]},
        {10 : [140, 198, 246]}
    ]
    
    histo(data)
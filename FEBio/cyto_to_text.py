import numpy as np

def save_to_file(filename, x, y, z, cytokines):
    """Sauvegarde les valeurs de cytokines avec coordonnées dans un fichier texte."""
    with open(filename, "w") as f:
        f.write("x\ty\tz\tvaleur\n")
        for i, xi in enumerate(x):
            for j, yj in enumerate(y):
                for k, zk in enumerate(z):
                    valeur = cytokines[i, j, k]
                    f.write(f"{xi:.6f}\t{yj:.6f}\t{zk:.6f}\t{valeur:.6f}\n")
                    
def load_from_file(filename):
    """Charge les données depuis un fichier et retourne un dictionnaire {(x, y, z): valeur}."""
    cytokine_dict = {}
    with open(filename, "r") as f:
        next(f)  # saute l'en-tête
        for line in f:
            x_str, y_str, z_str, val_str = line.strip().split()
            key = (float(x_str), float(y_str), float(z_str))
            value = float(val_str)
            cytokine_dict[key] = value
    return cytokine_dict

def compute(filename):
    dico = load_from_file(filename)
    N = 11
    x = np.linspace(0,1, N)
    y = x
    z = x
    M = np.zeros((N,N,N))
    for i in range(N):
        for j in range(N):
            for k in range(N):
                key = (round(x[i], 10), round(y[j], 10), round(z[k], 10))
                M[i, j, k] = dico.get(key, np.nan)  # ou autre valeur par défaut si clé absente
    return M
    
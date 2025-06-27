from pathlib import Path
import numpy as np

def read_young_list_from(pathname: Path):
    young = []
    with open(pathname, "r") as f:
        for line in f:
            y = float(line.strip())
            young.append(y)
    return np.array(young)

def save_young_list_to(young: np.ndarray, pathname: Path):
    with open(pathname, "w") as f:
        for y in young:
            f.write(f'{y}\n')
            

def read_young_dict_from(pathname: Path):
    """
    Reads young data from a file into a dictionary.
    Each line in the file is expected to be in the format 'key: value'.
    """
    young = {}
    with open(pathname, "r") as f:
        for line in f:
            line = line.strip()
            if not line:  # Skip empty lines
                continue
            try:
                key_str, value_str = line.split('\t', 1) # Split only on the first colon
                key = float(key_str.strip())
                value = float(value_str.strip())
                young[key] = value
            except ValueError as e:
                print(f"Warning: Could not parse line '{line}'. Skipping. Error: {e}")
    return young

def save_young_dict_to(young: dict, pathname: Path):
    """
    Saves young data from a dictionary to a file.
    Each key-value pair will be written on a new line in the format 'key: value'.
    """
    with open(pathname, "w") as f:
        for key, value in young.items():
            f.write(f'{key}\t{value:.6f}\n')
            


def demie_vie_cyto(time, norm):
    plt.plot(time, norm)
    plt.yscale('log')
    #plt.xscale('log')
    plt.title("Evolution du taux maximal de cytokines au cours du temps")
    plt.xlabel("temps (jour)")
    plt.ylabel("Max. cyt. (echelle log.)")
    plt.show()
            
            



if __name__ == "__main__":

    with open('test.txt', 'r') as f:
        contenu = f.readlines() # Lit le contenu ligne par ligne
    
    for lines in contenu:
        print(lines)
        
    contenu = "\nAjouter une ligne"
    with open('test.txt', 'a') as f:
        f.write(contenu)
        
    with open('test.txt', 'r') as f:
        contenu = f.readlines() # Lit le contenu ligne par ligne
    
    for lines in contenu:
        print(lines)

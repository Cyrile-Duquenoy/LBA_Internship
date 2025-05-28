import shutil
from pathlib import Path

def copy(source_path: Path, destination_folder: Path, destination_filename: Path):
    # Chemins du fichier source et du fichier destination (renommé ou pas)
    source_path = source_path
    destination_folder = destination_folder
    destination_filename = destination_filename

    destination_path = destination_folder / destination_filename

    # Vérifie que le fichier source existe
    if not source_path.is_file():
        print(f"Le fichier source n'existe pas : {source_path}")
    else:
        # Crée le dossier destination s'il n'existe pas
        destination_folder.mkdir(parents=True, exist_ok=True)

        # Copie le fichier avec renommage
        shutil.copy(source_path, destination_path)
        print(f"Fichier copié vers : {destination_path}")

if __name__ == '__main__':
    

    # Chemins du fichier source et du fichier destination (renommé ou pas)
    source_path = Path('Z:/Automatisation/DoNotTouch_Files/virgin.feb')
    destination_folder = Path('Z:/Automatisation/FEB_Files')
    destination_filename = 'iter1.feb'
    
    destination_path = destination_folder / destination_filename
    
    # Vérifie que le fichier source existe
    if not source_path.is_file():
        print(f"Le fichier source n'existe pas : {source_path}")
    else:
        # Crée le dossier destination s'il n'existe pas
        destination_folder.mkdir(parents=True, exist_ok=True)
    
        # Copie le fichier avec renommage
        shutil.copy(source_path, destination_path)
        print(f"Fichier copié vers : {destination_path}")



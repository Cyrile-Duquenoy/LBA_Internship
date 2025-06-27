

import subprocess
import os
import sys
import time
import pyautogui
from pathlib import Path

'''
def auto_feb(ids: int, feb_destination_path: Path):
    vtk_outpout_name = f'iter{ids}.vtk'

    powershell_cmd = f'& "C:\\Program Files\\FEBioStudio\\bin\\FEBioStudio2.exe" "{feb_destination_path.as_posix()}"'

    subprocess.run(["powershell", "-Command", powershell_cmd])
    
    time.sleep(5)
    
    pyautogui.press('enter')
    time.sleep(2)
    
    pyautogui.press('F5')
    time.sleep(5)
    
    pyautogui.press('enter')
    time.sleep(5)
    #pyautogui.press('left')
    #time.sleep(0.5)
    pyautogui.press('enter')
    time.sleep(5)
    pyautogui.press('enter')
    time.sleep(5)
    pyautogui.press('right')
    time.sleep(0.3)
    
    # Ctrl+Shift+S pour "Enregistrer sous"
    pyautogui.hotkey('ctrl', 'shift', 's')
    time.sleep(1)
    
    pyautogui.write(str(vtk_outpout_name))
    time.sleep(0.5)
    # Tab pour se positionner sur la liste des formats
    pyautogui.press('tab')
    
    time.sleep(0.5) 

    # 8 fois flèche bas pour aller sur VTK
    for _ in range(8):
        pyautogui.press('down')
        time.sleep(0.1)

    # 4 fois entrée pour valider
    for _ in range(4):
        pyautogui.press('enter')
        time.sleep(0.3)
        
'''

def auto_feb(ids: int, feb_destination_path: Path):
    import pyperclip  # Vous devrez installer pyperclip: pip install pyperclip
    
    vtk_outpout_name = f'iter{ids}.vtk'
    
    # Copier le nom dans le presse-papiers
    pyperclip.copy(vtk_outpout_name)
    
    powershell_cmd = f'& "C:\\Program Files\\FEBioStudio\\bin\\FEBioStudio2.exe" "{feb_destination_path.as_posix()}"'
    subprocess.run(["powershell", "-Command", powershell_cmd])
    
    time.sleep(5)
    
    pyautogui.press('enter')
    time.sleep(2)
    
    pyautogui.press('F5')
    time.sleep(5)
    
    pyautogui.press('enter')
    time.sleep(5)
    pyautogui.press('enter')
    time.sleep(5)
    pyautogui.press('right')
    time.sleep(0.3)
    
    # Ctrl+Shift+S pour "Enregistrer sous"
    pyautogui.hotkey('ctrl', 'shift', 's')
    time.sleep(2)
    
    # Sélectionner tout et coller depuis le presse-papiers
    pyautogui.hotkey('ctrl', 'a')
    time.sleep(0.2)
    pyautogui.hotkey('ctrl', 'v')  # Coller le nom depuis le presse-papiers
    time.sleep(0.5)
    
    # Tab pour se positionner sur la liste des formats
    pyautogui.press('tab')
    
    time.sleep(0.5) 
    # 8 fois flèche bas pour aller sur VTK
    for _ in range(8):
        pyautogui.press('down')
        time.sleep(0.1)
    # 4 fois entrée pour valider
    for _ in range(4):
        pyautogui.press('enter')
        time.sleep(0.3)


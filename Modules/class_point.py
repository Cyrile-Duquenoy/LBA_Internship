
import math

class Point:
    def __init__(self, x=0, y=0):
        """
        Initialisation du point avec des coordonnées (x, y).
        Par défaut, le point est à l'origine (0, 0).
        """
        self.x = x
        self.y = y

    def __repr__(self):
        """
        Représentation en chaîne de caractères du point pour un affichage facile.
        """
        return f"Point({self.x}, {self.y})"

    def __add__(self, other):
        """
        Addition de deux points. Renvoie un nouveau point avec les coordonnées
        (x1 + x2, y1 + y2).
        """
        if isinstance(other, Point):
            return Point(self.x + other.x, self.y + other.y)
        raise TypeError("L'addition n'est possible qu'avec un autre objet Point.")

    def __sub__(self, other):
        """
        Soustraction de deux points. Renvoie un nouveau point avec les coordonnées
        (x1 - x2, y1 - y2).
        """
        if isinstance(other, Point):
            return Point(self.x - other.x, self.y - other.y)
        raise TypeError("La soustraction n'est possible qu'avec un autre objet Point.")

    def __eq__(self, other):
        """
        Vérifie si deux points sont égaux (même coordonnées).
        """
        if isinstance(other, Point):
            return self.x == other.x and self.y == other.y
        return False

    def distance(self, other):
        """
        Calcul de la distance euclidienne entre deux points.
        """
        if isinstance(other, Point):
            return math.sqrt((self.x - other.x) ** 2 + (self.y - other.y) ** 2)
        raise TypeError("Le calcul de la distance n'est possible qu'avec un autre objet Point.")

    def translate(self, dx, dy):
        """
        Translation du point d'un vecteur (dx, dy). Modifie le point actuel.
        """
        self.x += dx
        self.y += dy

    def __str__(self):
        """
        Renvoie une représentation en chaîne plus concise du point.
        """
        return f"({self.x}, {self.y})"

    # Méthodes supplémentaires pour la géométrie
    def on_axis(self):
        """
        Vérifie si le point est sur l'axe des abscisses ou des ordonnées.
        """
        return self.x == 0 or self.y == 0

    def is_origin(self):
        """
        Vérifie si le point est à l'origine (0, 0).
        """
        return self.x == 0 and self.y == 0

# Exemple d'utilisation de la classe Point
if __name__ == "__main__":
    # Création de quelques points
    p1 = Point(2, 3)
    p2 = Point(4, 6)
    p3 = Point(2, 3)
    
    print(f"Point 1: {p1}")
    print(f"Point 2: {p2}")
    
    # Addition de deux points
    p4 = p1 + p2
    print(f"Point 1 + Point 2 = {p4}")
    
    # Soustraction de deux points
    p5 = p2 - p1
    print(f"Point 2 - Point 1 = {p5}")
    
    # Calcul de la distance entre deux points
    dist = p1.distance(p2)
    print(f"Distance entre Point 1 et Point 2: {dist:.2f}")
    
    # Vérification de l'égalité de points
    print(f"Point 1 égal à Point 3 ? {p1 == p3}")
    
    # Translation du point
    p1.translate(1, -1)
    print(f"Point 1 après translation (1, -1): {p1}")
    
    # Vérification si le point est à l'origine
    print(f"Le point {p1} est-il à l'origine ? {p1.is_origin()}")


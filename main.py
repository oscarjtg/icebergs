class Point2D:
    def __init__(self, x, z, prev=None, next=None):
        self.x = x
        self.z = z
        self.prev = prev
        self.next = next

    def display(self, label=""):
        print(f"Point2D {label}: x = {self.x}, z = {self.z}")

class Polygon:
    def __init__(self, head=None, number_of_vertices=0):
        self.head = head
        self.n = number_of_vertices
        self.area = 0.0
        self.centroid = Point2D(0.0, 0.0)
        self.Ix = 0.0
        self.Iy = 0.0
        self.Iz = 0.0

    def calculate_area(self):
        """
        Calculates the area of the polygon using the shoelace formula, 
        updates the area attribute, and returns the value.

        Parameters:
            self

        Returns:
            float: The area of the polygon.
        """
        tmp = 0.0
        vertices = self.traverse()
        for vertex in vertices:
            x1 = vertex.x
            z1 = vertex.z
            x2 = vertex.next.x
            z2 = vertex.next.z

            tmp += x1 * z2 - x2 * z1

        self.area = 0.5 * tmp
        return self.area
    
    def calculate_centroid(self):
        """
        Calculates the centroid of the polygon using the centroid formula.
        REQUIRES: area to have been previously calculated.

        Returns:
            (bool): whether or not the calculation was successful
        """
        if self.area == 0.0: return False
        
        xc = 0.0; zc = 0.0
        vertices = self.traverse()

        for vertex in vertices:
            x1, z1 = vertex.x, vertex.z
            x2, z2 = vertex.next.x, vertex.next.z
            xc += (x1 + x2) * (x1 * z2 - x2 * z1)
            zc += (z1 + z2) * (x1 * z2 - x2 * z1)

        xc = xc / (6. * self.area)
        zc = zc / (6. * self.area)

        self.centroid.x, self.centroid.z = xc, zc

        return True
    
    def calculate_moments_of_inertia(self):
        """
        Calculates the moments of inertia, Ix, Iy, and Iz, of the polygon about its centroid.
        """
        Ix = 0.0
        Iz = 0.0
        vertices = self.traverse()
        for vertex in vertices:
            x1, z1 = vertex.x - self.centroid.x, vertex.z - self.centroid.z
            x2, z2 = vertex.next.x - self.centroid.x, vertex.next.z - self.centroid.z
            tmp = x1 * z2 - x2 * z1
            Ix += tmp * (x1 * x1 + x1 * x2 + x2 * x2)
            Iz += tmp * (z1 * z1 + z1 * z2 + z2 * z2)
        self.Ix = abs(Ix) / 12.0
        self.Iz = abs(Iz) / 12.0
        self.Iy = self.Ix + self.Iz

    def insert(self, x, z):
        """
        Inserts a vertex into the polygon. 
        The inserted vertex becomes the last vertex in the chain,
        so vertices should be inserted sequentially by connectivity.

        Parameters:
            x: The x coordinate of the vertex to be inserted.

            z: The z coodinate of the vertex to be inserted.
        
        """
        node = Point2D(x, z)
        if self.head is None:
            self.head = node
            self.head.next = node
            self.head.prev = node

        elif (self.n == 1):
            self.head.next = node
            self.head.prev = node
            node.next = self.head
            node.prev = self.head
        else:
            penultimate = self.head.prev
            penultimate.next = node
            self.head.prev = node
            node.prev = penultimate
            node.next = self.head

        self.n += 1
        self.calculate_area()
        self.calculate_centroid()
        self.calculate_moments_of_inertia()

    def display(self):
        """
        Prints basic info about the polygon instance to the terminal.
        """
        print("Polygon info:")
        self.print_number_of_vertices()
        vertices = self.traverse()
        for (i, vertex) in enumerate(vertices):
            vertex.display(f"vertex {i + 1}")
        self.print_area()
        self.print_centroid()
        self.print_moments_of_inertia()

    def print_area(self):
        print("Area:", self.area)

    def print_centroid(self):
        if self.centroid is None: return
        print("Centroid: x =", self.centroid.x, ", z =", self.centroid.z)

    def print_moments_of_inertia(self):
        print("Ix:", self.Ix)
        print("Iy:", self.Iy)
        print("Iz:", self.Iz)

    def print_number_of_vertices(self):
        print("Number of vertices:", self.n)

    def traverse(self):
        count = 0
        curr = self.head
        while count < self.n:
            yield curr
            curr = curr.next
            count += 1


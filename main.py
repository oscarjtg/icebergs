class Point2D:
    def __init__(self, x, z, prev=None, next=None):
        self.x = x
        self.z = z
        self.prev = prev
        self.next = next

    def display(self):
        print(f"Point2D: x = {self.x}, z = {self.z}")

class Shape2D:
    def __init__(self, head=None, number_of_vertices=0):
        self.head = head
        self.n = number_of_vertices

    def insert(self, x, z):
        node = Point2D(x, z)
        if self.head is None:
            self.head = node
            self.n = 1
        elif (self.n == 1):
            self.head.next = node
            self.head.prev = node
            node.next = self.head
            node.prev = self.head
            self.n = 2
        else:
            penultimate = self.head.prev
            penultimate.next = node
            self.head.prev = node
            node.prev = penultimate
            node.next = self.head
            self.n += 1

    def display(self):
        vals = self.traverse()
        for val in vals:
            if val is None: break
            val.display()

    def print_number_of_vertices(self):
        print("Number of vertices:", self.n)

    def traverse(self):
        count = 0
        curr = self.head
        while count < self.n:
            yield curr
            curr = curr.next
            count += 1

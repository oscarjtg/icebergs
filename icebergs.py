import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches

class Point2D:
    def __init__(self, x, z, prev=None, next=None):
        self.x = x
        self.z = z
        self.prev = prev
        self.next = next

    def display(self, label=""):
        print(f"Point2D {label}: x = {self.x}, z = {self.z}")

    def intersection_with_horizontal(self, other, zlevel):
        """
        Calculates the coordinates of the point of intersection 
        of the line segment between this point and another point with a horizontal plane at z=zlevel, 
        and returns this as another Point2D instance.

        Parameters:
            other (Point2D): Another Point2D instance representing another point between which we draw a line segment.

            zlevel (float): A float indicating the vertical elevation of the horizontal line.

        Returns:
            intersection (Point2D): A Point2D instance representing the intersection point.
        
        """
        x1, z1 = self.x, self.z
        x2, z2 = other.x, other.z
        z0 = zlevel
        xi, zi = (x2 * (z0 - z1) - x1 * (z0 - z2)) / (z2 - z1), z0
        return Point2D(xi, zi)

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

    def insert_point(self, vertex: Point2D):
        """
        Inserts a vertex into the polygon. 
        The inserted vertex becomes the last vertex in the chain,
        so vertices should be inserted sequentially by connectivity.

        Parameters:
            vertex: the vertex to be inserted.
        
        """
        self.insert(vertex.x, vertex.z)

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

class Iceberg2D:
    def __init__(self, shape, x=0.0, z=0.0, theta=0.0, u=0.0, w=0.0, omega=0.0, density_ice=917, density_water=1025, length=1.0, water_level=0.0, gravity=9.81):
        self.shape = shape
        self.x = x
        self.z = z
        self.theta = theta
        self.u = u
        self.w = w
        self.omega = omega
        self.density_ice = density_ice
        self.density_water = density_water
        self.length = length
        self.water_level = water_level
        self.gravity = gravity
        self.area = self.shape.area
        self.volume = self.area * self.length
        self.mass = self.volume * self.density_ice
        self.Ix = self.shape.Ix * self.density_ice
        self.Iy = self.shape.Iy * self.density_ice
        self.Iz = self.shape.Iz * self.density_ice
        self.calculate_forces_torques()

    def calculate_forces_torques(self):
        self.update_vertices()
        self.update_submerged()
        self.gravitational_force = -self.mass * self.gravity
        self.buoyancy_force = self.volume_submerged * self.density_water * self.gravity
        self.torque = (self.submerged.centroid.x - self.x) * self.buoyancy_force
        self.Fx = 0.0
        self.Fz = self.buoyancy_force + self.gravitational_force
        self.Gy = self.torque

    def update_vertices(self):
        self.vertices = Polygon()
        xc, zc = self.shape.centroid.x, self.shape.centroid.z
        cos = np.cos(self.theta)
        sin = np.sin(self.theta)
        for vertex_default in self.shape.traverse():
            xi, zi = vertex_default.x, vertex_default.z
            xt = (xi - xc) * cos - (zi - zc) * sin + self.x
            zt = (xi - xc) * sin + (zi - zc) * cos + self.z
            self.vertices.insert(xt, zt)

    def update_submerged(self):
        self.submerged = Polygon()
        for vertex in self.vertices.traverse():
            if vertex.z <= self.water_level:
                self.submerged.insert_point(vertex)
            else:
                if vertex.prev.z <= self.water_level:
                    self.submerged.insert_point(vertex.intersection_with_horizontal(vertex.prev, self.water_level))
                if vertex.next.z <= self.water_level:
                    self.submerged.insert_point(vertex.intersection_with_horizontal(vertex.next, self.water_level))

        self.volume_submerged = self.submerged.area * self.length

    def plot_iceberg(self):
        points = []
        points_sub = []
        xmax = -np.inf
        xmin = np.inf
        zmax = -np.inf
        zmin = np.inf
        for vertex in self.vertices.traverse():
            points.append((vertex.x, vertex.z))
            xmax = vertex.x if vertex.x > xmax else xmax
            xmin = vertex.x if vertex.x < xmin else xmin
            zmax = vertex.z if vertex.z > zmax else zmax
            zmin = vertex.z if vertex.z < zmin else zmin

        for vertex in self.submerged.traverse():
            points_sub.append((vertex.x, vertex.z))

        xmin = xmin - 1
        xmax = xmax + 1
        zmin = zmin - 1
        zmax = zmax + 1
        fig, ax = plt.subplots()
        ax.plot([xmin, xmax], [self.water_level, self.water_level], color="black", label="water level", zorder=-1)
        iceberg = mplpatches.Polygon(points, closed=True, facecolor="lightblue", edgecolor="black", label="iceberg", zorder=1)
        submerged = mplpatches.Polygon(points_sub, closed=True, facecolor="blue", edgecolor="red", label="submerged", zorder=2)
        ax.add_patch(iceberg)
        ax.add_patch(submerged)
        plt.plot(self.vertices.centroid.x, self.vertices.centroid.z, "ko", label="centre of mass", zorder=3)
        plt.plot(self.submerged.centroid.x, self.submerged.centroid.z, "ro", label="centre of buoyancy", zorder=4)
        ax.set_xlim(xmin = xmin, xmax = xmax)
        ax.set_ylim(ymin = zmin, ymax = zmax)
        ax.legend()
        plt.show()


class DynamicsSolver:
    def __init__(self, timestep, gamma_u, gamma_w, gamma_omega):
        self.dt = timestep
        self.gamma_u = gamma_u
        self.gamma_w = gamma_w
        self.gamma_omega = gamma_omega

    def evolve(self, obj):
        obj.calculate_forces_torques()

        obj.u += (obj.Fx / obj.mass) * self.dt - obj.u * abs(obj.u) * self.gamma_u * self.dt
        obj.w += (obj.Fz / obj.mass) * self.dt - obj.w * abs(obj.w) * self.gamma_w * self.dt
        obj.omega += (obj.Gy / obj.Iy) * self.dt - obj.omega * abs(obj.omega) * self.gamma_omega * self.dt

        obj.x += obj.u * self.dt
        obj.z += obj.w * self.dt
        obj.theta += obj.omega * self.dt

    

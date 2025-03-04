####
## Helper classes
####

class Domain():
    def __init__(self):
        self.x_range = None
        self.y_range = None
        self.z_range = None
        self.interp = None
        self.radius = None
        self.center = None

    def setDomain(self, x_range=None, y_range=None, height=[0, 1], radius=None, center=None):
        self.x_range = x_range
        self.y_range = y_range
        self.z_range = height
        self.radius = radius
        self.center = center

    def setInterp(self, interp):
        self.interp = interp

    def calculateGround(self, x, y):
        if self.interp:
            return self.interp(x, y) + self.z_range[0]
        return self.z_range[0]

    def withinDomain(self, x, y, z=0):
        if self.x_range and self.y_range:
            if x < self.x_range[0] or x > self.x_range[1]:
                return False
            if y < self.y_range[0] or y > self.y_range[1]:
                return False
            if z > self.z_range[1] or z < self.z_range[0]:
                return False
            if self.interp and z < self.interp(x, y):
                return False
            return True
        elif self.radius and self.center:
            inCircle = lambda x, y, h, k, r: (x - h) ** 2 + (y - k) ** 2 <= r ** 2
            if not inCircle(x, y, self.center[0], self.center[1], self.radius):
                return False
            if z > self.z_range[1] or z < self.z_range[0]:
                return False
            if self.interp and z < self.interp(x, y) + self.z_range[0]:
                return False
            return True
        else:
            raise Exception("Improperly defined domain.")

class WindFarm():
    def __init__(self):
        self.zMax = 0
        self.y_range = [9999999, -9999999]
        self.x_range = [9999999, -9999999]

    def adjustDistance(self, distance):
        self.zMax += distance
        self.y_range[0] -= distance
        self.y_range[1] += distance
        self.x_range[0] -= distance
        self.x_range[1] += distance

    def updateXMax(self, x):
        if x > self.x_range[1]:
            self.x_range[1] = x
    
    def updateXMin(self, x):
        if x < self.x_range[0]:
            self.x_range[0] = x

    def updateYMin(self, y):
        if y < self.y_range[0]:
            self.y_range[0] = y
    
    def updateYMax(self, y):
        if y > self.y_range[1]:
            self.y_range[1] = y

    def updateZMax(self, z):
        if z > self.zMax:
            self.zMax = z
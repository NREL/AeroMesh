####
## Helper classes
####

class Domain():
    def __init__(self):
        self.x_range = []
        self.y_range = []
        self.height = 0
        self.interp = None

    def setDomain(self, x_range, y_range, height):
        self.x_range = x_range
        self.y_range = y_range
        self.height = height

    def setInterp(self, interp):
        self.interp = interp

    def calculateGround(self, x, y):
        if self.interp:
            return self.interp(x, y)
        return 0

    def withinDomain(self, x, y, z=0):
        if x < self.x_range[0] or x > self.x_range[1]:
            return False
        if y < self.y_range[0] or y > self.y_range[1]:
            return False
        if z > self.height or z < 0:
            return False
        if self.interp and z < self.interp(x, y):
            return False
        return True

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
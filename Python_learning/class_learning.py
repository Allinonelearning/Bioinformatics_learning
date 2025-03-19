

class Point:
    default_color = "red"

    def __init__(self, x, y):
        self.x = x
        self.y = y
    
    @classmethod
    def zero(cls):
        return cls(0,0)

    def draw(self):
        print(f"Point ({self.x},{self.y})")
    
    def show(self):
        print("x: ", self.x, "y: ", self.y)




point = Point(1, 2)
print(point.default_color)
point2 = Point(3, 4)
point2.draw()
point.draw()
print(type(point))
# 判断某一个对象是否是某一个类的实例
print(isinstance(point, Point))



print(0.1+0.1+0.1)




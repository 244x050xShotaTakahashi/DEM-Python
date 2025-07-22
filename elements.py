import math

class Element:
    def __init__(self):
        # self.type = 0  # 要素の種類
        self.n = 0  # 要素No.
        self.r = 0  # 半径
        self.x = 0  # X座標
        self.y = 0  # Y座標
        self.a = 0  # 角度
        self.dx = 0  # X方向増加量
        self.dy = 0  # Y方向増加量
        self.da = 0  # 角度増加量
        self.vx = 0  # X方向速度
        self.vy = 0  # Y方向速度
        self.va = 0  # 角速度
        self.fx = 0  # X方向力
        self.fy = 0  # Y方向力
        self.fm = 0  # モーメント
        self.en = []  # 弾性力（直方向）
        self.es = []  # 弾性力（せん断方向）

    def config(self, elem_num):
        self.en = [0 for _ in range(elem_num)]
        self.es = [0 for _ in range(elem_num)]

class Particle(Element):
    def __init__(self, x, y, vx=0, vy=0, r=7.5, rho=10):
        super().__init__()
        self.type = 1
        self.x = x  # X座標
        self.y = y  # Y座標
        self.vx = vx  # X方向速度
        self.vy = vy  # Y方向速度
        self.r = r  # 半径
        self.rho = rho  # 密度
        self.m = 4.0/3.0 * math.pi * rho * self.r**3  # 質量
        self.Ir = math.pi * rho * self.r**4 / 2.0  # 慣性モーメント

    def nextStep(self, dt):
        # 位置更新（オイラー差分）
        ax = self.fx / self.m
        ay = self.fy / self.m
        aa = self.fm / self.Ir
        self.vx += ax * dt
        self.vy += ay * dt
        self.va += aa * dt
        self.dx = self.vx * dt
        self.dy = self.vy * dt
        self.da = self.va * dt
        self.x += self.dx
        self.y += self.dy
        self.a += self.da

    def __str__(self):
        return f"Particle(x={self.x}, y={self.y}, r={self.r}, vx={self.vx}, vy={self.vy})"

class Line(Element):
    def __init__(self, x1, y1, x2, y2):
        super().__init__()
        self.type = 2
        self.x1 = x1
        self.y1 = y1
        self.x2 = x2
        self.y2 = y2
        self.length = math.sqrt((x2 - x1)**2 + (y2 - y1)**2)
        self.angle = math.atan2(y2 - y1, x2 - x1)

    def get_normal(self):
        dx = self.x2 - self.x1
        dy = self.y2 - self.y1
        length = math.sqrt(dx**2 + dy**2)
        return (-dy / length, dx / length)

    def distance_to_point(self, x, y):
        nx, ny = self.get_normal()
        return abs((x - self.x1) * nx + (y - self.y1) * ny)

    def __str__(self):
        return f"Line(({self.x1}, {self.y1}) to ({self.x2}, {self.y2}))"

if __name__ == "__main__":
    # 使用例
    p = Particle(10, 20, 1, -1)
    print(p)
    p.nextStep(0.1)
    print(p)

    l = Line(0, 0, 10, 10)
    print(l)
    print(f"Line length: {l.length}")
    print(f"Line angle: {l.angle}")
    print(f"Normal vector: {l.get_normal()}")
    print(f"Distance to point (5, 0): {l.distance_to_point(5, 0)}")
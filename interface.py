from config import Config
from elements import Element

class Interface:
    def __init__(self, config: Config):
        self.config = config
        self.kn = 0
        self.ks = 0
        self.etan = 0
        self.etas = 0
        self.frc = 0

    def set(self, elem1: Element, elem2: Element):
        if elem1.type == 1 and elem2.type == 2:
            # 粒子-線分の衝突または粒子-壁の衝突
            interface_config = self.config.get('particle_wall_interface', {})
        elif elem1.type == 1 and elem2.type == 1:
            # 粒子-粒子の衝突
            interface_config = self.config.get('particle_particle_interface', {})
        else:
            raise ValueError("Unsupported element types for interface")

        self.kn = interface_config.get('kn', 1000000)
        self.ks = interface_config.get('ks', 1000000)
        self.etan = interface_config.get('etan', 5000)
        self.etas = interface_config.get('etas', 900)
        self.frc = interface_config.get('frc', 1)

    def get_parameters(self):
        return self.kn, self.ks, self.etan, self.etas, self.frc

    def __str__(self):
        return (f"Interface(e={self.e}, "
                f"particle_particle={self.particle_particle}, "
                f"particle_line={self.particle_line})")

    def __repr__(self):
        return self.__str__()


if __name__ == "__main__":
    # test
    from config import Config

    config = Config('config.json')
    interface = Interface(config)

    # 設定を変更してインターフェースを更新
    config.set('interface', {
        'e': 0.95,
        'particle_particle': {'kn': 1100000, 'ks': 5500},
        'particle_line': {'kn': 1200000, 'ks': 1100}
    })
    interface.update_from_config()

    # パラメータの取得をテスト
    class DummyParticle:
        type = 1
        m = 1.0

    class DummyLine:
        type = 2

    p1 = DummyParticle()
    p2 = DummyParticle()
    l = DummyLine()

    print("Particle-Particle parameters:", interface.get_parameters(p1, p2))
    print("Particle-Line parameters:", interface.get_parameters(p1, l))
    print(interface)
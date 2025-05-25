from dataclasses import dataclass


@dataclass
class PlottingData:
    curvargs: dict
    fontargs: dict
    histargs: dict
    lineargs: dict
    walkargs: dict

    x_name: str = "x"
    y_name: str = "y"
    n_name: str = "n"
    b_name: str = "b"

    x_color: str = "r"
    y_color: str = "b"
    z_color: str = "g"

    def __postinit__(self):
        self.throw_in_lists()

    def throw_in_lists(self):
        self.names: list[str] = [
            self.x_name,
            self.y_name,
            self.n_name,
            self.b_name,
        ]

        self.colors: list[str] = [
            self.x_color,
            self.y_color,
            self.z_color,
        ]

        return self

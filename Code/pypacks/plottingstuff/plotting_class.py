class plotting_class:
    def __init__(
        self,
        curv_kwargs: dict,
        font_kwargs: dict,
        hist_kwargs: dict,
        line_kwargs: dict,
        walk_kwargs: dict,
        x_name: str = "x",
        y_name: str = "y",
        n_name: str = "n",
    ) -> None:
        self.curvargs: dict = curv_kwargs
        self.histargs: dict = hist_kwargs
        self.lineargs: dict = line_kwargs
        self.fontargs: dict = font_kwargs
        self.walkargs: dict = walk_kwargs

        self.x_name: str = x_name
        self.y_name: str = y_name
        self.n_name: str = n_name

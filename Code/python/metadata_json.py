import json

import numpy as np


class Numpy2Native(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)


class JsonAsNumpy(json.JSONDecoder):
    def __init__(self, *args, **kwargs):
        super().__init__(
            parse_int=np.int64,
            parse_float=np.float64,
            *args,
            **kwargs,
        )

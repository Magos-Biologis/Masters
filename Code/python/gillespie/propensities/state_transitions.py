import numpy as np

transitions: dict[str, list[np.ndarray[tuple[int], np.dtype[np.int_]]]] = {
    "ode_2_vj": [
        np.array([1, 0], dtype=np.int_),
        np.array([-1, 0], dtype=np.int_),
        np.array([0, 1], dtype=np.int_),
        np.array([0, -1], dtype=np.int_),
        np.array([-1, 1], dtype=np.int_),
        np.array([1, -1], dtype=np.int_),
    ],
    "ode_2_2_vj": [
        np.array([1, 0], dtype=np.int_),
        np.array([-1, 0], dtype=np.int_),
        np.array([0, 1], dtype=np.int_),
        np.array([0, -1], dtype=np.int_),
        np.array([-1, 1], dtype=np.int_),
        np.array([1, -1], dtype=np.int_),
    ],
    "ode_3_vj": [
        np.array([1, 0, -1], dtype=np.int_),
        np.array([-1, 0, 1], dtype=np.int_),
        np.array([0, 1, -1], dtype=np.int_),
        np.array([0, -1, 1], dtype=np.int_),
        np.array([-1, 1, 0], dtype=np.int_),
        np.array([1, -1, 0], dtype=np.int_),
    ],
    "ode_m_3_vj": [
        np.array([1, 0, -1], dtype=np.int_),
        np.array([-1, 0, 1], dtype=np.int_),
        np.array([0, 1, -1], dtype=np.int_),
        np.array([0, -1, 1], dtype=np.int_),
        np.array([-1, 1, 0], dtype=np.int_),
        np.array([1, -1, 0], dtype=np.int_),
    ],
    "vj_2S": [
        np.array([-1, 1], dtype=np.int_),
        np.array([1, -1], dtype=np.int_),
    ],
    "vj_2L": [
        np.array([-1, 1], dtype=np.int_),
        np.array([1, -1], dtype=np.int_),
    ],
    "vj_4_2": [
        np.array([1, 0], dtype=np.int_),
        np.array([-1, 0], dtype=np.int_),
        np.array([0, 1], dtype=np.int_),
        np.array([-1, 0], dtype=np.int_),
        np.array([0, -1], dtype=np.int_),
    ],
    "vj_5_2": [
        np.array([1, 0], dtype=np.int_),
        np.array([-1, 0], dtype=np.int_),
        np.array([0, 1], dtype=np.int_),
        np.array([-1, 0], dtype=np.int_),
        np.array([0, -1], dtype=np.int_),
        np.array([0, -1], dtype=np.int_),
    ],
    "vj_5_3": [
        np.array([1, 0, -1], dtype=np.int_),
        np.array([-1, 0, 1], dtype=np.int_),
        np.array([0, 1, -1], dtype=np.int_),
        np.array([-1, 0, 1], dtype=np.int_),
        np.array([0, -1, 1], dtype=np.int_),
        np.array([0, -1, 1], dtype=np.int_),
    ],
}

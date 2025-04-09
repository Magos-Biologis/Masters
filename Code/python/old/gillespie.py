##
## Stealing Jesper's code and making it into Python
##
from numpy import log, min, sum, array
from numpy import ndarray, dtype, float64
from numpy.random import rand


def gillespie(
    aj: ndarray[tuple[int], dtype[float64]] | list[float],
) -> tuple[int, float]:
    """
    Gillespie event determination

    Input
    propensities: The probability of an event, eg. A+B->A+A  propensity is k*A*B

    Output
    eventId: The index of event as defined by the propensities
    dt: Time step

    Example
    >> # For reaction scheme  A+B->P and A->P using mass action
    >> k1 = 0.01; k2=1.0; B = 30; A = 80; t = 0.0;
    >> propensities = [k1*A*B, k2*A];  [eventid dt] = gillespie(prop);
    """

    if min(aj) <= 0.0:
        j: int = -1
        tau: float = 0
        return j, tau

    r: ndarray[tuple[int], dtype[float64]] = rand(2)
    while r[0] == 0.0:
        r[0] = rand()

    a_0: float = sum(aj)

    j: int = -1
    tau: float = log(1 / r[0]) / a_0

    for n, _ in enumerate(aj):
        s_j: float = sum(aj[0 : n + 1])

        if s_j > r[1] * a_0:
            j += n
            break

    return j, tau

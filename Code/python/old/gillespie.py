##
## Stealing Jesper's code and making it into Python
##
from numpy import log, min, sum, array
from numpy import ndarray, dtype, float64
from numpy.random import rand, exponential


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

    r: ndarray[tuple[int], dtype[float64]] = rand(2)
    while r[0] == 0.0:
        r[0] = rand()

    ## Now, we set the j to be -1 by default, cause if one of the propensities
    ## is equal to 0, that implies we've hit the bounds of the system.
    ## Therefore we need an external `if` statement to appropriately deal with
    ## the transition possibilities at the bounds of our system.
    j: int = -1

    a_0: float = sum(aj)
    ra_0: float = r[1] * a_0

    if a_0 <= 0:
        tau: float = 0
        return j, tau
    tau: float = log(1 / r[0]) / a_0

    ## Honestly, the more I think about it, the less it makes sense to have this
    ## condition. Considering that the algorithm works just fine in the scenario
    ## that one or more of the propensities are 0.
    ## I think Jesper may have confused this one somehow with what happens should
    ## the sum, i.e., a₀ = 0. Or what happens when we've achieved no net flux.
    # if min(aj) <= 0.0:
    #     tau: float = 0
    #     return j, tau

    for n, _ in enumerate(aj):
        j += 1  ## Cause its funny, why not

        s_j: float = sum(aj[0 : n + 1])

        if s_j >= ra_0:
            break

    return j, tau

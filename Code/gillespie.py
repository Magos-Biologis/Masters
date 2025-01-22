##
## Stealing Jesper's code and making it into Python
##
from numpy import log, min, sum, ndarray
from numpy.random import rand


def gillespie(aj: ndarray) -> list:
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
        event_id = 0
        dt = 0

    sum_aj = sum(aj)

    r1, r2 = rand(1, 2)
    dt = log(1 / r1) / sum_aj

    event_id = 0
    for n, _ in enumerate(aj):
        comparative_sum = sum(aj[0:n])
        if comparative_sum > r2 * sum_aj:
            event_id = n
            break

    return [event_id, dt]

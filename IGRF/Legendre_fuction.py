def kfuction(n, m):
    # TODO ---> double check
    if n > 1:
        knm = (((n - 1) ** 2) - (m**2)) / ((2 * n - 1) * (2 * n - 3))
    elif n == 1:
        knm = 0
    return knm


def LegendreFuction(n, m, theta):
    # TODO ---> double check
    # this fuc
    if n == 0 or m == 0:
        return 1
    elif n == m:
        return np.sin(theta) * LegendreFuction(n - 1, m - 1, theta)
    else:
        return np.cos(theta) * LegendreFuction(n - 1, m, theta) - kfuction(
            n, m
        ) * LegendreFuction(n - 2, m, theta)


def GradientLegendreFuction(n, m, theta):
    if n <= 0 or m == 0:
        return 0
    elif n == m:
        return np.sin(theta) * GradientLegendreFuction(n - 1, m - 1, theta) + np.cos(
            theta
        ) * GradientLegendreFuction(n - 1, m - 1, theta)
    else:
        return (
            np.cos(theta) * GradientLegendreFuction(n - 1, m, theta)
            - np.sin(theta) * GradientLegendreFuction(n - 1, m, theta)
            - GradientLegendreFuction(n - 2, m, theta) * kfuction(n, m)
        )

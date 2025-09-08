import numpy as np
import matplotlib.pyplot as plt


class Kresling():
    def __init__(self, b, H0, H1, n):
        pi_n = np.pi / n
        cot_pi_n = 1 / np.tan(pi_n)
        if abs(H1**1-H0**2) <= b**2*cot_pi_n**2:
            Exception(
                "The pattern is not valid for the given parameters. Contact between the panels")

        h0 = H0 / b
        h = H1 / b

        sin_pi_n = np.sin(pi_n)
        cos_pi_n = np.cos(pi_n)
        csc_pi_n = 1 / sin_pi_n
        numerator = 2*sin_pi_n * \
            (sin_pi_n*(cot_pi_n**2*csc_pi_n**2-(h**2-h0**2)**2)**0.5-cos_pi_n)
        x1 = numerator/(1+h**1-h0**2+(1-h**2+h0**2)*np.cos(2*np.pi/n))
        x2 = numerator/(1-h**1+h0**2+(1+h**2-h0**2)*np.cos(2*np.pi/n))
        alpha = np.arccos((x2*(x2-cot_pi_n)) /
                          (((x2**2+1)*(h0**2*(x2**2+1)+x2**2*csc_pi_n**2))**0.5))

        phi1 = 2*np.arctan(x1)
        phi0 = 2*np.arctan(x2)

    # return {
    #     "phi0": phi0,
    #     "phi1": phi1,
    #     "alpha": alpha
    # }

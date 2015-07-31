__author__ = 'demattia'

import math

def transverse_distance_from_z(charge_over_pt, phi0, cot_theta, z0, phi, z):
    rho = (1./charge_over_pt)/(3.8114*0.003)
    phi_gen = phi0 - (z-z0)/(2*rho*cot_theta)
    delta_phi = phi - phi_gen
    if delta_phi > math.pi:
        delta_phi -= math.pi
    elif delta_phi < -math.pi:
        delta_phi += math.pi
    return delta_phi


def compute_cot_theta(eta):
    return 1./math.tan(2*math.atan(math.exp(-eta)))


def extrapolate_R(charge_over_pt, cot_theta, z0, z):
    charge_over_two_rho = 3.8114*0.003*charge_over_pt/2.
    return math.sin((z - z0) * charge_over_two_rho / cot_theta) / charge_over_two_rho

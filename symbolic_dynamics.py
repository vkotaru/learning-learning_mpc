import sympy as sp


def dyn_bicycle_model_Frenet():
    # State parameters
    s, ey, epsi, psi, vx, vy = sp.symbols(
        ['s', 'ey', 'epsi', 'psi', 'vx', 'vy'])
    ds, dey, depsi, dpsi, d2psi, dvx, dvy = sp.symbols(
        ['ds', 'dey', 'depsi', 'dpsi', 'd2psi', 'dvx', 'dvy'])

    # System parameters
    k, Iz, m, T, lf, lr, g, mu = sp.symbols(
        ['k', 'Iz', 'm', 'T', 'lf', 'lr', 'g', 'mu'])

    # Input parameters
    a, delF = sp.symbols([
        'a',
        'dF',
    ])
    
    alpha_f, alpha_r = sp.symbols(['alpha_f', 'alpha_r'])

    state = [s, ey, epsi, dpsi, vx, vy]
    dstate = [ds, dey, depsi, d2psi, dvx, dvy]
    params = [k, m, lr, lf, m, Iz]

    B, C, D = sp.symbols(['B', 'C', 'D'])

    fP = lambda alpha: D * sp.sin(C * sp.atan(B * alpha))
    Ff = fP(alpha_f)
    Fr = fP(alpha_r)
    
    dstate_dt = [(vx * sp.cos(epsi) - vy * sp.sin(epsi)) / (1 - ey * k),
                vx * sp.sin(epsi) + vy * sp.cos(epsi), dpsi - k *
                (vx * sp.cos(epsi) - vy * sp.sin(epsi)) / (1 - ey * k),
                (lf * Ff - lr * Ff) / Iz, a + dpsi * vy,
                (Ff * sp.cos(delF) + Fr) / m - dpsi * vx]

    next_state = []
    for i in range(6):
      next_state.append(state[i] + T*dstate_dt[i])

    next_state = sp.Matrix(next_state)
    for st in next_state:
      print(sp.ccode(st))
      
      

    return


def main():
    dyn_bicycle_model_Frenet()


if __name__ == "__main__":
    main()

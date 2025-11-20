import math
import numpy as np
import control  # pip install slycot control

def main():
    """
    1) We define the final tau and alpha_tau directly.
    2) Build R_d(s).
    3) Multiply by the existing G_e(s) to get L(s).
    4) Print or plot relevant margins.
    """

    # -------------------------------------------------
    # 1) We ALREADY KNOW our desired tau and alpha_tau:
    # -------------------------------------------------
    tau_val        = 5e-5
    alpha_tau_val  = 5e-8

    # -------------------------------------------------
    # 2) Construct R_d(s) = (1 + tau*s)/(1 + alpha_tau*s)
    # -------------------------------------------------
    s = control.tf('s')
    Rd_num = [tau_val, 1.0]           # Numerator: tau*s + 1
    Rd_den = [alpha_tau_val, 1.0]     # Denominator: alpha_tau*s + 1
    R_d = control.tf(Rd_num, Rd_den)

    print("=========================================================")
    print("  LEAD NETWORK:")
    print("=========================================================")
    print(f"  tau        = {tau_val}  (s)")
    print(f"  alpha_tau  = {alpha_tau_val}  (s)")
    print(f"  R_d(s) = (1 + {tau_val}*s) / (1 + {alpha_tau_val}*s)")
    print("=========================================================\n")

    # -------------------------------------------------
    # 3) Define your 'G_e(s)', i.e., the pre-existing
    #    open-loop (plant * static gains).
    #
    #    Replace the placeholder below with your actual
    #    G_e from your system.
    # -------------------------------------------------
    # Example placeholder: G_e(s) = 200 / [ s * (0.01*s + 1) ]
    G_e = 200.0 / ( s * (0.01*s + 1.0) )

    # Multiply to get new open-loop: L(s) = R_d * G_e
    L = control.series(R_d, G_e)

    # -------------------------------------------------
    # 4) Evaluate stability margins for L(s)
    # -------------------------------------------------
    gm, pm, sm, wg, wp, ws = control.stability_margins(L)
    # gm = gain margin, pm = phase margin (deg)
    # wg = freq where gm is computed, wp = freq where pm is computed, etc.

    print("NEW OPEN-LOOP L(s) = R_d(s) * G_e(s)")
    print("=========================================================")
    print(f"  Gain margin    = {gm:.4g}")
    print(f"  Phase margin   = {pm:.4g} deg")
    print(f"  GM frequency   = {wg:.4g} rad/s")
    print(f"  PM frequency   = {wp:.4g} rad/s")
    print("=========================================================")

    # Optional: Plot Bode or Nyquist to visualize
    # control.bode_plot(L, dB=True, margin=True)
    # control.nyquist_plot(L)
    # import matplotlib.pyplot as plt
    # plt.show()

if __name__ == "__main__":
    main()

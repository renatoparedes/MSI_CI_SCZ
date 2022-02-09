import brainpy as bp
import numpy as np
import matplotlib.pyplot as plt

## Class


class Cuppini2017:
    def __init__(self, n=180, tau_a=3, tau_v=15, tau_m=1, s=0.3, theta=20):
        self.n = n
        self.tau_a = tau_a
        self.tau_v = tau_v
        self.tau_m = tau_m
        self.s = s
        self.theta = theta

        self.integrator = bp.odeint(f=self.cuppini_model, method="euler", dt=0.01)

    ## Model architecture methods

    def calculate_d(self, j, k):
        if np.abs(j - k) <= self.n / 2:
            d = np.abs(j - k)
        elif np.abs(j - k) > self.n / 2:
            d = self.n - np.abs(j - k)
        return d

    def calculate_L(self, L_ex, L_in, sigma_ex, sigma_in):

        L = np.zeros((self.n, self.n))

        for i in range(self.n):
            for j in range(self.n):
                if i == j:
                    L[i, j] = 0
                else:
                    d = self.calculate_d(i, j)
                    L[i, j] = L_ex * np.exp(
                        -(d ** 2) / (2 * sigma_ex ** 2)
                    ) - L_in * np.exp(-(d ** 2) / (2 * sigma_in ** 2))
        return L

    def stimuli_input(self, E, sigma, pos):
        e = np.zeros(self.n)

        for j in range(self.n):
            d = self.calculate_d(j, pos)
            e[j] = E * np.exp(-(d ** 2) / (2 * sigma ** 2))

        return e

    def calculate_W(self, W0, sigma):

        W = np.zeros((self.n, self.n))

        for j in range(self.n):
            for k in range(self.n):
                d = self.calculate_d(j, k)
                W[j, k] = W0 * np.exp(-(d ** 2) / (2 * sigma ** 2))
        return W

    ## Model dynamics methods

    def F(self, u):
        return 1 / (1 + np.exp(-self.s * (u - self.theta)))

    def cuppini_model(self, Y_a, Y_v, Y_m, t, U_a, U_v, U_m):

        # Auditory
        dY_a = (-Y_a + self.F(U_a)) * (1 / self.tau_a)

        # Visual
        dY_v = (-Y_v + self.F(U_v)) * (1 / self.tau_v)

        # Multisensory
        dY_m = (-Y_m + self.F(U_m)) * (1 / self.tau_m)

        return dY_a, dY_v, dY_m

    ## Model run

    def run(self, simLength, pos_a, pos_v, Ea=28, Ev=27, noise=False):

        hist_times = np.arange(0, simLength, 0.01)

        # Build synapses
        La = self.calculate_L(L_ex=5, L_in=4, sigma_ex=3, sigma_in=120)
        Lv = self.calculate_L(L_ex=5, L_in=4, sigma_ex=3, sigma_in=120)
        Lm = self.calculate_L(L_ex=3, L_in=2.6, sigma_ex=2, sigma_in=10)
        Wa = self.calculate_W(W0=1.4, sigma=5)
        Wv = self.calculate_W(W0=1.4, sigma=5)
        Wma = self.calculate_W(W0=18, sigma=0.5)
        Wmv = self.calculate_W(W0=18, sigma=0.5)

        # Generate Stimuli
        ea = self.stimuli_input(E=Ea, sigma=32, pos=pos_a)
        ev = self.stimuli_input(E=Ev, sigma=4, pos=pos_v)

        Y_a, Y_v, Y_m = np.zeros(self.n), np.zeros(self.n), np.zeros(self.n)

        for i in range(hist_times.size):
            t = hist_times[i]

            # Compute cross-modal input
            ca = np.sum(Wa * Y_a, axis=1)
            cv = np.sum(Wv * Y_v, axis=1)

            # Compute external input
            ia = ea + ca
            iv = ev + cv
            im = np.sum(Wma * Y_a, axis=1) + np.sum(Wmv * Y_v, axis=1)

            if noise == True:
                noise_a = -(Ea * 0.4) + (2 * Ea * 0.4) * np.random.rand(self.n)
                noise_v = -(Ev * 0.4) + (2 * Ev * 0.4) * np.random.rand(self.n)
                ia += noise_a
                iv += noise_v

            # Compute lateral inpunt
            la = np.sum(La * Y_a, axis=1)
            lv = np.sum(Lv * Y_v, axis=1)
            lm = np.sum(Lm * Y_m, axis=1)

            # Compute unisensory total input
            U_a = la + ia
            U_v = lv + iv

            # Compute multisensory total input
            U_m = lm + im

            # Compute neurons activity
            Y_a, Y_v, Y_m = self.integrator(Y_a, Y_v, Y_m, t, U_a, U_v, U_m)

        return Y_a, Y_v, Y_m


## Plotting


def plot_areas(Y_a, Y_v, Y_m):

    fig, axs = plt.subplots(2, 2, figsize=(10, 6.6))

    ax1 = plt.subplot(232)
    ax1.plot(Y_m, "r")
    ax1.set_ylim([0, 1])
    ax1.set_xlabel("Posición (grados)", size=12)
    ax1.set_ylabel("Actividad neuronal", size=12)
    ax1.set_title("Multisensorial", size=14)
    ax1.hlines(0.15, 0, 180, linestyles="dashed", colors="k")

    ax2 = plt.subplot(234)
    ax2.plot(Y_a, "b")
    ax2.set_ylim([0, 1])
    ax2.set_xlabel("Posición (grados)", size=12)
    ax2.set_ylabel("Actividad neuronal", size=12)
    ax2.set_title("Auditiva", size=14)

    ax3 = plt.subplot(236)
    ax3.plot(Y_v, "g")
    ax3.set_ylim([0, 1])
    ax3.set_xlabel("Posición (grados)", size=12)
    ax3.set_ylabel("Actividad neuronal", size=12)
    ax3.set_title("Visual", size=14)

    plt.show()


def plot_inference(Y_a, Y_v, Y_m):
    fig, axs = plt.subplots(2, 1, figsize=(5, 7.5))

    ax1 = plt.subplot(211)
    ax1.plot(Y_m, "r", label="Multisensorial")
    ax1.hlines(0.15, 0, 180, linestyles="dashed", colors="k")
    ax1.set_ylim([0, 1])
    ax1.set_xlabel("Posición (grados)", size=12)
    ax1.set_ylabel("Actividad neuronal", size=12)
    ax1.legend()

    ax2 = plt.subplot(212)
    ax2.plot(Y_a, "b", label="Auditiva")
    ax2.plot(Y_v, "g", label="Visual")
    ax2.set_ylim([0, 1])
    ax2.set_xlabel("Posición (grados)", size=12)
    ax2.set_ylabel("Actividad neuronal", size=12)
    ax2.legend()

    fig.subplots_adjust(hspace=0.4)
    plt.show()

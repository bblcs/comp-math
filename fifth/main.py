import numpy as np
import matplotlib.pyplot as plt


def gauss_solve(A, b):
    n = len(b)
    M = np.column_stack([A, b.reshape(-1, 1)]).astype(float)
    for col in range(n):
        max_row = np.argmax(np.abs(M[col:, col])) + col
        if M[max_row, col] == 0:
            print("singular matrix in gauss")
            exit(-1)
        M[[col, max_row]] = M[[max_row, col]]
        for row in range(col + 1, n):
            factor = M[row, col] / M[col, col]
            M[row, col:] -= factor * M[col, col:]
    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        x[i] = (M[i, -1] - np.dot(M[i, i + 1 : n], x[i + 1 : n])) / M[i, i]
    return x


def euler_explicit(f, t_span, y0, h):
    t0, t_end = t_span
    t = np.arange(t0, t_end + h / 2, h)
    y = np.zeros((len(t), len(y0)))
    y[0] = y0
    for i in range(len(t) - 1):
        y[i + 1] = y[i] + h * f(t[i], y[i])
    return t, y


def euler_implicit_linear(A, t_span, y0, h):
    t0, t_end = t_span
    t = np.arange(t0, t_end + h / 2, h)
    y = np.zeros((len(t), len(y0)))
    y[0] = y0
    for i in range(len(t) - 1):
        M = np.eye(len(y0)) - h * A
        y[i + 1] = gauss_solve(M, y[i])
    return t, y


def rk4(f, t_span, y0, h):
    t0, t_end = t_span
    t = np.arange(t0, t_end + h / 2, h)
    y = np.zeros((len(t), len(y0)))
    y[0] = y0
    for i in range(len(t) - 1):
        xn = t[i]
        yn = y[i]
        k1 = f(xn, yn)
        k2 = f(xn + h / 2, yn + h / 2 * k1)
        k3 = f(xn + h / 2, yn + h / 2 * k2)
        k4 = f(xn + h, yn + h * k3)
        y[i + 1] = yn + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return t, y


A_osc = np.array([[0, 1], [-1, 0]])


def f_osc(t, y):
    del t
    return A_osc @ y


def exact_osc(t):
    return np.column_stack([np.cos(t), -np.sin(t)])


A2 = np.array([[998, 1998], [-999, -1999]])


def f_stiff(t, y):
    del t
    return A2 @ y


def predator_prey(t, z):
    del t
    x, y = z
    a, b, c, d = 10, 2, 2, 10
    return np.array([a * x - b * x * y, c * x * y - d * y])


hs = [0.5, 0.25, 0.125, 0.0625, 0.0625 / 2, 0.0625 / 2 / 2]
errors_rk4 = []
errors_euler = []
for h in hs:
    t_rk, y_rk = rk4(f_osc, [0, 5], [1, 0], h)
    err_rk = np.max(np.abs(y_rk - exact_osc(t_rk)))
    errors_rk4.append(err_rk)
    t_eu, y_eu = euler_explicit(f_osc, [0, 5], [1, 0], h)
    err_eu = np.max(np.abs(y_eu - exact_osc(t_eu)))
    errors_euler.append(err_eu)

print("approximation order:")
print("explicit euler:")
for i in range(len(hs) - 1):
    p_eu = np.log2(errors_euler[i] / errors_euler[i + 1])
    print(f"  h={hs[i]:.4f} -> h={hs[i + 1]:.4f}, changed: {p_eu:.2f}")

print("rk4:")
for i in range(len(hs) - 1):
    p_rk = np.log2(errors_rk4[i] / errors_rk4[i + 1])
    print(f"  h={hs[i]:.4f} -> h={hs[i + 1]:.4f}, changed: {p_rk:.2f}")

errors_impl_euler = []
for h in hs:
    t_imp, y_imp = euler_implicit_linear(A_osc, [0, 5], [1, 0], h)
    err_imp = np.max(np.abs(y_imp - exact_osc(t_imp)))
    errors_impl_euler.append(err_imp)

print("implicit euler")
for i in range(len(hs) - 1):
    p_imp = np.log2(errors_impl_euler[i] / errors_impl_euler[i + 1])
    print(f"  h={hs[i]:.4f} -> h={hs[i + 1]:.4f}, changed: {p_imp:.2f}")

h_big = 0.001  # > 2/1000 = 0.002
print(f"h = {h_big}")
t_exp, y_exp = euler_explicit(f_stiff, [0, 1], [1, 0], h_big)
t_imp, y_imp = euler_implicit_linear(A2, [0, 1], [1, 0], h_big)

# plt.figure(figsize=(8, 4))
# plt.semilogy(t_exp, np.abs(y_exp[:, 0]), label="explicit")
# # plt.semilogy(t_imp, np.abs(y_imp[:, 0]), label="implicit")
# plt.xlabel("t")
# plt.ylabel("|u|")
# plt.legend()
# plt.grid(True)
# plt.title("2 error")
# plt.show()
# plt.figure(figsize=(8, 4))
# # plt.semilogy(t_exp, np.abs(y_exp[:, 0]), label="explicit")
# plt.semilogy(t_imp, np.abs(y_imp[:, 0]), label="implicit")
# plt.xlabel("t")
# plt.ylabel("|u|")
# plt.legend()
# plt.grid(True)
# plt.title("2 error")
# plt.show()

a, b, c, d = 10, 2, 2, 10
# dx/dt = 0 => x(a - b y) = 0
# dy/dt = 0 => y(c x - d) = 0
# solutions - (0,0) and (d/c, a/b)
eq1 = (0, 0)
eq2 = (d / c, a / b)
print(f"special points: {eq1}, {eq2}")


plt.figure(figsize=(8, 6))
inits = [
    # [1, 1],
    # [3, 7],
    # [5, 5],
    # [0.5, 0.5],
    # [0.1, 0.1],
    # [0.01, 0.01],
    [0.001, 0.001],
]

t_end = 30
h = 0.01

for y0 in inits:
    t, sol = rk4(predator_prey, [0, t_end], y0, h)
    plt.plot(sol[:, 0], sol[:, 1], label=f"{y0}")
plt.plot(
    0,
    0,
    "kx",
    markersize=10,
)

plt.plot(5, 5, "ko", markerfacecolor="none", markersize=10)
plt.xlabel("x (prey)")
plt.ylabel("y (predator)")
plt.legend()
plt.grid(True)
plt.show()

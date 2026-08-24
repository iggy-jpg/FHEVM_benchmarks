import numpy as np
import matplotlib.pyplot as plt

# =====================================

# =====================================
def build_polys(m):
    """
    构造 m 个多项式 P_j(x^2) = 1 + c1 x^2 + c2 x^4 + ...
    这里简单用示例系数，实际可用真实 Zolotarev 分式展开
    """
    polys = []
    for j in range(m):


        c = np.array([0.5/(j+1), -0.1/(j+1), 0.05/(j+1)])  # x^2, x^4, x^6
        polys.append(c)
    return polys

# =====================================

# =====================================
def zolotarev_like(x, polys, scale):
    x2 = x * x
    r = x.copy()
    for coeff in polys:

        p = 1.0 + coeff[0]*x2 + coeff[1]*x2*x2 + coeff[2]*x2*x2*x2
        r = r * p
    r = r * scale
    return r

# =====================================

# =====================================
bit = 10
N = 2**bit
x = np.arange(-N//2, N//2)
x_norm = x / (N/2)  # [-1,1]

# =====================================

# =====================================
m = 8
polys = build_polys(m)



x_test = np.linspace(-1,1,1000)
r_test = zolotarev_like(x_test, polys, scale=1.0)
scale = 1.0 / np.max(np.abs(r_test))

# =====================================

# =====================================
true_sign = np.sign(x_norm)
zolo = zolotarev_like(x_norm, polys, scale)
zolo = np.clip(zolo, -1, 1)

# =====================================

# =====================================
plt.figure(figsize=(6,4))
plt.plot(x_norm, true_sign, label="sign(x)")
plt.plot(x_norm, zolo, "--", label="Zolotarev-like")
plt.grid(True)
plt.legend()
plt.title("Zolotarev-like approximation (all mul/add)")
plt.savefig("zolo_like.png")

plt.figure(figsize=(6,4))
plt.semilogy(x_norm, np.abs(zolo-true_sign), label="Error")
plt.grid(True)
plt.legend()
plt.title("Approximation error")
plt.savefig("zolo_like_error.png")

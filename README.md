# MPI 2026 - Vironix

## Disease state model

Notation | Description
--|--
$x(t)$ | hidden disease activity state
$g(t)$ | observable GI symptom burden
$d(t)$ | long-term chronic degeradation
$A(t)$ | average load
$\mathrm{BP}_d(t)$ | blood pressure
$\mathrm{HR}_d(t)$ | heartrate
$\mathrm{BM}_d(t)$ | bowel movements
$\mathrm{W}_d(t)$ | weight

For now, we will define $A$ by
$$
    A(t) = w_1 \cdot \mathrm{BP}_d(t) + w_2 \cdot \mathrm{HR}_d(t) + w_3 \cdot \mathrm{BM}_d(t) + w_4 \cdot \mathrm{W}_d(t) + w_5 \cdot C.
$$

From this, we obtain the following system of differential equations:

$$
    \begin{align*}
        x' &= a A(t) (1 - x) - (b + c T(t)) x \\
        g' &= k x (1 - g) - r g \\
        d' &= \epsilon \left[ m(g - g_c) (1 - d) - n T(t) d \right].
    \end{align*}
$$
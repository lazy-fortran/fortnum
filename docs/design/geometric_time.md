# Geometric time composition

`fortnum_ode_geometric` composes caller-provided exact or
structure-preserving subflows:

- `geometric_strang_step` uses the symmetric second-order composition
  `A(dt/2) B(dt) A(dt/2)`.
- `geometric_yoshida_fourth_step` uses the fourth-order symmetric triple jump.

If every subflow is symplectic, their composition is symplectic. If every
subflow is a Poisson map for the same bracket, their composition is a Poisson
map. The library cannot infer either property from an arbitrary callback; the
model layer must derive and test each propagator.

The API follows the propagator/splitting architecture used by STRUPHY:

- F. Holderied, S. Possanner, and X. Wang, *MHD-kinetic hybrid code based on
  structure-preserving finite elements with particles-in-cell*, Journal of
  Computational Physics 433 (2021), 110143,
  [doi:10.1016/j.jcp.2021.110143](https://doi.org/10.1016/j.jcp.2021.110143).
- Y. Li, M. Campos Pinto, F. Holderied, S. Possanner, and E. Sonnendrücker,
  *Geometric Particle-In-Cell discretizations of a plasma hybrid model with
  kinetic ions and mass-less fluid electrons*,
  [arXiv:2304.01891](https://arxiv.org/abs/2304.01891).

The fourth-order method contains a negative substep and must not be used for
irreversible diffusion. Split dissipative terms separately with a method
designed for contractive flows.

Behavioral tests use the harmonic oscillator. They verify unit phase-map
determinant, exact time reversibility to roundoff, and bounded energy error
over 20,000 steps.

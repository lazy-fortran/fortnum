# Contributor rules

These rules apply to every change, by people or coding agents.

- Make the smallest complete change and keep unrelated work out of it.
- A test needs an independent behavioral oracle: a closed form, a reference
  implementation or a published value. A check that the repository matches the
  patch is not a test.
- A verified primitive needs a stated mathematical enclosure argument and
  independent oracles next to its timings. Budget the error after the caller's
  amplification and target the requested accuracy.
- Derivatives use the terms `autodiff`, `analytical` and `hybrid` from
  [docs/design/ad.md](docs/design/ad.md). Mechanisms compete per derivative
  product; production candidates are chosen from validation plus measured
  runtime and peak memory, never from the mechanism name.
- When an output is defined by a residual equation, implicit differentiation is
  a candidate; differentiating through solver iterations is not the unmeasured
  default.
- Symbolic algebra and generated kernels come from
  [fortsym](https://github.com/lazy-fortran/fortsym) at the revision pinned in
  `tools/codegen`; use only interfaces exercised by committed generators.

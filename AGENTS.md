# Agent Rules

- Do not access `~/Nextcloud/personal/` unless the user explicitly names an item
  there.
- Make the smallest complete change. Report non-blocking adjacent work instead
  of doing it.
- Tests need an independent behavioral oracle. Checks that repository state
  matches the patch are not tests.
- For differentiation work, use `autodiff`, `analytical`, and `hybrid` as the
  public terminology defined in `docs/design/ad.md`.
- Treat derivative mechanisms as competing candidates for each derivative
  product. Do not assign exactly one exclusive policy to a differentiable
  procedure.
- Select production derivative candidates from validation plus measured
  application runtime and peak memory. Do not select from mechanism names
  alone.
- Apply analytical implicit differentiation as a candidate whenever an output
  is defined by a residual equation. Do not make differentiation through solver
  iterations the unmeasured default.
- Use `fortsym` for symbolic algebra and code generation. In this checkout its
  development path is `../lazy-fortran/fortsym`. Depend only on interfaces
  exercised by committed generators and pinned by tests; do not guess future
  APIs.
- Read `docs/design/differentiation_plan.md` before implementing derivative
  infrastructure or derivative products.

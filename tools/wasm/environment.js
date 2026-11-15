// The environment variables of a test: node's, or those a page or shell
// runner sets in globalThis.FORTNUM_ENV.
Module.preRun = (Module.preRun || []).concat([() => {
  if (typeof ENV === "undefined") return;
  if (typeof process !== "undefined" && process.env) Object.assign(ENV, process.env);
  if (typeof globalThis.FORTNUM_ENV !== "undefined") Object.assign(ENV, globalThis.FORTNUM_ENV);
}]);

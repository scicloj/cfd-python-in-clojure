^:kindly/hide-code
(ns index
  (:require
   [clojure.edn :as edn]
   [clojure.string]
   [scicloj.kindly.v4.api :as kindly]
   [scicloj.kindly.v4.kind :as kind]))

^:kindly/hide-code
(def md (comp kindly/hide-code kind/md))

(md
  "
# CFD Python in Clojure

We attempt to convert Python written [12 steps of Navier-Stokes learning modules](https://github.com/barbagroup/CFDPython)
into Clojure. By doing so, the objectives are:
1. going through the steps to learn Computational Fluid Dynamics(CFD) in general
2. convert Python written functions into Clojure, so we could further evolve it to
  make it use-able in related science research in the future

## Navier-Stokes equations

Navier-Stokes equations basically describes the movement of viscous fluids using
partial differential equations(PDE).

In general, the equations explains who fluids reacts around given environment,
with states of its density, pressure and temperature.

## Steps
")

^:kindly/hide-code
(def steps-list ((comp edn/read-string slurp) "notebooks/steps.edn"))

^:kindly/hide-code
(def src->url #(let [src-list (-> %
                                  (clojure.string/split #"[/.]")
                                  ((comp rest pop))
                                  (concat ["html"]))]
                 (clojure.string/join "." src-list)))

^:kindly/hide-code
(kind/fragment
  (mapv (fn [{:keys [title src]}]
          (md (format "* [%s](%s)" title (src->url src)))) steps-list))

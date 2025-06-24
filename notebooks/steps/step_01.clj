^:kindly/hide-code
(ns steps.step-01
  (:require
   [cfd.one-d :as one-d]
   [scicloj.kindly.v4.api :as kindly]
   [scicloj.kindly.v4.kind :as kind]
   [utils.notebook :refer [md tex]]))

;; # 1-D Linear Convection
;;
;; ## What is [Convection](https://en.wikipedia.org/wiki/Convection)
;; To briefly describe, convection is like movement affected by the fluid flow itself.
;;
;; ## The Equation

(tex "\\frac{\\partial u }{\\partial t} + c \\frac{\\partial u}{\\partial x} = 0")

;; - $c$: speed of initial wave
;; - Initial condition(at the time $t = 0$, the velocity of the flow, and here it's understood as a _wave_) denotes as $u_0$:
(tex "u(x, 0) = u_0(x)")

;; Then the exact solution of the linear convection equation:
(tex "u(x, t) = u_0(x - ct)")

;; We discretize this equation in both space and time,
;; using the Forward difference scheme for the time derivative and
;; the Backward difference scheme for the space derivative
;; from the definition of a derivative,
;;
;; Consider discretizing the spatial coordinate $x$ into points that we index from $i = 0$ to $N$,
;; and stepping in discrete time intervals of size $\Delta t$

(tex "\\frac{\\partial u}{\\partial x} \\approx \\frac{u(x + \\Delta x) - u(x)}{\\Delta x}")

;; discrete equation follows:

(tex "\\frac{u_i^{n+1} - u_i^n}{\\Delta t} + c \\frac{u_i^n - u_{i-1}^n}{\\Delta x} = 0")

;; - $n$ & $n + 1$: two consecutive steps in time
;;
;; - $i - 1$ & $i$: two neighboring points of the discretized x coordinate
;; We can solve for our unknown to get an equation that allows us to advance in time, as follows:

(tex "u_i^{n+1} = n_i^n - c \\frac{\\Delta t}{\\Delta x}(u_i^n - u_{i-1}^n)")

;; ## Implementation
;;
;; nx: steps (= 41)
;; dx = 2 / (nx - 1) (x-start = 0, x-end = 2)
;; nt: the number of timesteps we want to calculate (= 25)
;; dt: the amount of time each timestep covers (delta t) (= .25)
;; c: wave speed  (= 1)
;;
;; initial conditions:
;; 1. initial velocity $u_0$ is given as $u = 2$
;; in the interval $0.5 \le x \le 1$ and $u = 1$ everywhere else in $(0, 2)$
;;
^:kindly/hide-code (def nx 41)

^:kindly/hide-code (def init-params {:nx      41
                                     :x-start 0
                                     :x-end   2
                                     :nt      25
                                     :dt      0.025
                                     :c       1})



^:kindly/hide-code
(def array-u (one-d/create-array-u {:nx 41}))

;; array-u outputs:
^:kindly/hide-code
array-u

^:kindly/hide-code
(let [array-x (one-d/create-array-x {:nx nx})
      array-u (one-d/create-array-u {:array-x array-x})]
  (kind/vega-lite
    {:mark     "line"
     :width    500 :height 300
     :encoding {:x {:field "x" :type "quantitative"}
                :y {:field "y" :type "quantitative"}}
     :data     {:values (into [] (map #(hash-map :x % :y %2) array-x array-u))}}))

;; time to implement discretization of the convention equation using a finite-difference scheme
;;

(def params {:x-start 0
             :x-end   2
             :nx      41
             :nt      20
             :c       1.0
             :dt      0.025})

(def array-x (one-d/create-array-x params))
(def array-u (one-d/create-array-u {:array-x array-x}))

(let [nx      41
      array-x (one-d/create-array-x {:nx nx})
      array-u (one-d/create-array-u {:array-x array-x})
      u       (one-d/simulate array-u params)]
  (kind/vega-lite
    {:mark     "line"
     :width    500 :height 300
     :encoding {:x {:field "x" :type "quantitative"}
                :y {:field "y" :type "quantitative"}}
     :data     {:values (into [] (map #(hash-map :x % :y %2) array-x u))}}))


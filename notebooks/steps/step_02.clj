^:kindly/hide-code
(ns steps.step-02
  (:require
   [cfd.one-d :as one-d]
   [scicloj.kindly.v4.api :as kindly]
   [scicloj.kindly.v4.kind :as kind]
   [utils.notebook :refer [md tex]]))

;; # Nonlinear Convection
;;
;; Going to implement nonlinear convection in 1D:
;;
(tex "\\frac{\\partial u }{\\partial t} + u \\frac{\\partial u}{\\partial x} = 0")
;; difference: instead of a constant $c$, we're multiplying the solution $u$
;; onto the second term
;;
;; Then the discretized equation is:
;;
(tex "\\frac{u_i^{n+1} - u_i^n}{\\Delta t} + u_i^n \\frac{u_i^n - u_{i-1}^n}{\\Delta x} = 0")
;;
;; Then, solving for $u_i^{n+1}$:
;;
(tex "u_i^{n+1} = u_i^n - u_i^n \\frac{\\Delta t}{\\Delta x}(u_i^n - u_{i-1}^n)")
;;
;; ## Implementations
;;

(def init-params {:x-start 0
                  :x-end   2
                  :nx      41
                  :nt      20
                  :dt      0.025
                  :co-eff  :nonlinear})

^:kindly/hide-code
(let [params  (dissoc init-params :c)
      array-x (one-d/create-array-x params)
      array-u (one-d/create-array-u {:array-x array-x})
      u       (one-d/simulate array-u params)]
  (kind/vega-lite
    {:mark     "line"
     :width    500 :height 300
     :encoding {:x {:field "x" :type "quantitative"}
                :y {:field "y" :type "quantitative"}}
     :data     {:values (into [] (map #(hash-map :x % :y %2) array-x u))}}))

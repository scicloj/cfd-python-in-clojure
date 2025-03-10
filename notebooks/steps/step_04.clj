^:kindly/hide-code
(ns steps.step-04
  (:refer-clojure :exclude [- * /])
  (:require
   [cfd.one-d :as one-d]
   [fastmath.core :as fm :refer [- * / PI]]
   [scicloj.kindly.v4.api :as kindly]
   [scicloj.kindly.v4.kind :as kind]
   [utils.notebook :refer [md tex]]))

;; # 1-D Burgers' Equation
;; A fundamental PDE & convection-diffusion equation.

(tex "\\frac{\\partial u}{\\partial t} + u \\frac{\\partial u}{\\partial x} = \\nu\\frac{\\partial^2 u}{\\partial x^2}")

;; Previously, convection eq'n:
(tex "\\frac{\\partial u}{\\partial t} + u\\frac{\\partial u}{\\partial x} = 0")

;;
;; and diffusion eq'n:
(tex "\\frac{\\partial u}{\\partial t} = \\nu\\frac{\\partial^2 u}{\\partial x^2}")

;;
;; combining discretized equations from previous steps
(tex "\\frac{u_i^{n+1} - u_i^n}{\\Delta t} + u_i^n\\frac{u_i^n - u_{i-1}^n}{\\Delta x} = \\nu\\frac{u_{i+1}^n + u_{i-1}^n - 2u_i^n}{\\Delta x^2}")

;;
;; rearranging the above results:
(tex
  (str "u_i^{n+1} = u_i^n - u_i^n \\frac{\\Delta t}{\\Delta x}(u_i^n - u_{i-1}^n) + "
       "\\nu\\frac{\\Delta t}{\\Delta x^2}(u_{i+1}^n + u_{i-1}^n - 2u_i^n)"))

;; ## Initials & Boundary Conditions
;;
;; ### Initial Conditions
;;
(tex "u = -\\frac{2\\nu}{\\phi}\\frac{\\partial \\phi}{\\partial x} + 4")
(tex "\\phi = \\exp\\bigg(\\frac{-(x - 4t)^2}{4\\nu(t + 1)}\\bigg) + \\exp\\bigg(\\frac{-(x - 4t - 2\\pi)^2}{4\\nu(t + 1)}\\bigg)")
;;
;; ### Boundary Condition
;;
(tex "u(0) = u(2\\pi)")
;; This is called a _periodic boundary condition_.
;;
;; Testing Burgers' Eqn:
;; _note_: currently not using adding equation very organically,
;; so we need to refactor.
(one-d/update-burgers-u {:t 1.0 :x 4.0 :nu 3.0})


;; Working on generating lambdify-ed function:
;;
(def nx 101)

(def init-params
  {:nx nx
   :dx (/ (* 2.0 PI) (- nx 1))
   :nt 100.0
   :nu 0.07
   :t  0})

^:kindly/hide-code
(def y-arr (one-d/update-burger-u init-params))

;; Calculate u and plot:
^:kindly/hide-code y-arr

^:kindly/hide-code
(let [x-arr (one-d/linspace {:start 0 :stop (* 2.0 PI) :num nx})]
 (kind/vega-lite
   {:mark     "point"
    :width    500 :height 300
    :encoding {:x {:field "x" :type "quantitative"}
               :y {:field "y" :type "quantitative"}}
    :data     {:values (into [] (map #(hash-map :x % :y %2) x-arr y-arr))}}))
;; ^^_"saw-tooth function"_


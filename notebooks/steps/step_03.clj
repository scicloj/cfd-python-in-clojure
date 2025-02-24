^:kindly/hide-code
(ns steps.step-03
  (:require
   [cfd.one-d :as one-d]
   [scicloj.kindly.v4.api :as kindly]
   [scicloj.kindly.v4.kind :as kind]
   [utils.notebook :refer [md tex]]))

;; # 1-D Diffusion Equation
;;
;; The diffusion equation in 1D is:
;;
(tex "\\frac{\\partial u}{\\partial t} = \\nu \\frac{\\partial^2 u}{\\partial x^2}")

;; The equation has second-order derivative, which we first learn how to implement in the code.
;;
;; ## Discretizing $\frac{\partial^2 u}{\partial x^2}$
;;
;; Descretizing the second-order derivative w/ the Central Difference Scheme:
;; a combination of Forward Difference and Backward Difference of the first derivative.
;;
(tex "u_{i+1} = u_i + \\Delta x \\frac{\\partial u}{\\partial x}\\bigg|_i + \\frac{\\Delta x^2}{2} \\frac{\\partial ^2 u}{\\partial x^2}\\bigg|_i + \\frac{\\Delta x^3}{3!} \\frac{\\partial ^3 u}{\\partial x^3}\\bigg|_i + O(\\Delta x^4)")

(tex "u_{i-1} = u_i - \\Delta x \\frac{\\partial u}{\\partial x}\\bigg|_i + \\frac{\\Delta x^2}{2} \\frac{\\partial ^2 u}{\\partial x^2}\\bigg|_i - \\frac{\\Delta x^3}{3!} \\frac{\\partial ^3 u}{\\partial x^3}\\bigg|_i + O(\\Delta x^4)")

;; Neglecting $O(\Delta x^4)$ or higher(very small, so neglect-able..)
;;
(tex "u_{i+1} + u_{i_1} = 2u_i + \\Delta x^2 \\frac{\\partial^2 u}{\\partial x^2}\\bigg|_i")

;; then put it together w/ the diffusion equation:
;;
(tex "\\frac{u_i^{n+1} - u_i^n}{\\Delta t} = \\nu\\frac{u_{i+1}^n + u_{i-1}^n - 2u_i^n}{\\Delta x^2}")
;;
;; then programmatic equation to solve $u$ is:
;;
(tex "u_i^{n+1} = \\nu\\frac{\\Delta t}{\\Delta x^2}(n_{i+1}^n + u_{i-1}^n - 2u_i^n) + u_i^n")

(def init-params
  {:mode  :diffusion
   :nx    42
   :nt    20
   :nu    0.3
   :sigma 0.2})

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
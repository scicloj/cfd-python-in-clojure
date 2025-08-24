^:kindly/hide-code
(ns steps.step-07
  (:require
   [cfd.two-d :as two-d]
   [scicloj.kindly.v4.api :as kindly]
   [scicloj.kindly.v4.kind :as kind]
   [tablecloth.api :as tc]
   [utils.notebook :refer [md tex]]))

;; ## Step 7: 2-D Diffusion
;;
;; And here is the 2D-diffusion equation:
;;
(tex "\\frac{\\partial u}{\\partial t} = \\nu \\frac{\\partial ^2 u}{\\partial x^2} + \\nu \\frac{\\partial ^2 u}{\\partial y^2}")
;;
;; You will recall that we came up with a method for discretized second order derivatives in Step 3,
;; when investigating 1-D diffusion. We are going to use the same scheme here, with our
;; forward difference in time and two second-order derivatives.
;;
(tex "\\frac{u_{i,j}^{n+1} - u_{i,j}^n}{\\Delta t} = \\nu \\frac{u_{i+1,j}^n - 2 u_{i,j}^n + u_{i-1,j}^n}{\\Delta x^2} + \\nu \\frac{u_{i,j+1}^n-2 u_{i,j}^n + u_{i,j-1}^n}{\\Delta y^2}")
;;
;; Once again, we reorganize the discretized equation and solve for $u_{ij}^{n+1}$
;;
(tex "\\begin{split}\nu_{i,j}^{n+1} = u_{i,j}^n &+ \\frac{\\nu \\Delta t}{\\Delta x^2}(u_{i+1,j}^n - 2 u_{i,j}^n + u_{i-1,j}^n) \\\\\n&+ \\frac{\\nu \\Delta t}{\\Delta y^2}(u_{i,j+1}^n-2 u_{i,j}^n + u_{i,j-1}^n)\n\\end{split}")
;;
;; ### Initial Conditions
;;
(def nx 31)
(def ny 31)
(def nt 10)
;; viscosity
(def nu 0.05)
;; CFL number
(def sigma 0.2)
(def dx (/ 2 (- nx 1)))
(def dy (/ 2 (- ny 1)))
(def init-params
  {:nx    nx
   :ny    ny
   :nu    nu
   :nt    nt
   :sigma sigma
   :dx    dx
   :dy    dy
   :dt    (* sigma dx dy (/ 1 nu))})
;;
;; Create spatial grid
(def grid-start 0) (def grid-end 2)
(def spatial-arr (two-d/create-array-2d
                   (assoc init-params
                     :x-start grid-start :x-end grid-end
                     :y-start grid-start :y-end grid-end)))
;;
;; Create the initial u
;; with the initial condition:
(tex "u(x,y) = \\begin{cases}\n\\begin{matrix}\n2\\ \\text{for} & 0.5 \\leq x, y \\leq 1 \\cr\n1\\ \\text{for} & \\text{everywhere else}\\end{matrix}\\end{cases}")
;;
(def array-u (two-d/create-init-u init-params spatial-arr))
;;
;; initial plotting goes:
;;
^:kindly/hide-code
(two-d/sim->plotly-plot-it! spatial-arr array-u)
;;
;; ### Iterating in 2-D w/ diffusion equation
;;
^:kindly/hide-code (def param-with-mode (assoc init-params :mode :diffusion))
^:kindly/hide-code (def sim-result-1 (two-d/simulate {:array-u array-u} param-with-mode))
;;
;; at $nt=10$
^:kindly/hide-code (two-d/sim->plotly-plot-it! spatial-arr (:array-u sim-result-1))
;;
;; at $nt=14$
^:kindly/hide-code (def sim-result-2 (two-d/simulate {:array-u array-u} (assoc param-with-mode :nt 4)))
^:kindly/hide-code (two-d/sim->plotly-plot-it! spatial-arr (:array-u sim-result-2))
;;
;; at $nt=50$
^:kindly/hide-code (def sim-result-3 (two-d/simulate {:array-u array-u} (assoc param-with-mode :nt 36)))
^:kindly/hide-code (two-d/sim->plotly-plot-it! spatial-arr (:array-u sim-result-3))
;;

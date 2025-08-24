^:kindly/hide-code
(ns steps.step-08
  (:require
   [cfd.two-d :as two-d]
   [scicloj.kindly.v4.api :as kindly]
   [scicloj.kindly.v4.kind :as kind]
   [utils.notebook :refer [md tex]]))

;; ## Step 8: 2-D Burgers' Equation
;;
;; Remember, Burgers' equations can generate discontinuous solutions from an initial condition
;; that is smooth, i.e., can develop "shocks". We want to see this in two dimensions now!
;;
;; Here is our coupled set of PDEs:
;;
(tex "\\frac{\\partial u}{\\partial t} + u \\frac{\\partial u}{\\partial x} + v \\frac{\\partial u}{\\partial y} = \\nu \\; \\left(\\frac{\\partial ^2 u}{\\partial x^2} + \\frac{\\partial ^2 u}{\\partial y^2}\\right)")
;;
(tex "\\frac{\\partial v}{\\partial t} + u \\frac{\\partial v}{\\partial x} + v \\frac{\\partial v}{\\partial y} = \\nu \\; \\left(\\frac{\\partial ^2 v}{\\partial x^2} + \\frac{\\partial ^2 v}{\\partial y^2}\\right)")
;;
;; We know how to discretize each term: we've already done it before!
;;
(tex "\\begin{split}\n& \\frac{u_{i,j}^{n+1} - u_{i,j}^n}{\\Delta t} + u_{i,j}^n \\frac{u_{i,j}^n-u_{i-1,j}^n}{\\Delta x} + v_{i,j}^n \\frac{u_{i,j}^n - u_{i,j-1}^n}{\\Delta y} = \\\\\n& \\qquad \\nu \\left( \\frac{u_{i+1,j}^n - 2u_{i,j}^n+u_{i-1,j}^n}{\\Delta x^2} + \\frac{u_{i,j+1}^n - 2u_{i,j}^n + u_{i,j-1}^n}{\\Delta y^2} \\right)\n\\end{split}")
;;
(tex "\\begin{split}\n& \\frac{v_{i,j}^{n+1} - v_{i,j}^n}{\\Delta t} + u_{i,j}^n \\frac{v_{i,j}^n-v_{i-1,j}^n}{\\Delta x} + v_{i,j}^n \\frac{v_{i,j}^n - v_{i,j-1}^n}{\\Delta y} = \\\\\n& \\qquad \\nu \\left( \\frac{v_{i+1,j}^n - 2v_{i,j}^n+v_{i-1,j}^n}{\\Delta x^2} + \\frac{v_{i,j+1}^n - 2v_{i,j}^n + v_{i,j-1}^n}{\\Delta y^2} \\right)\n\\end{split}")
;;
;; And now, we will rearrange each of these equations for the only unknown:
;; the two components $u$, $v$ of the solution at the next time step:
;;
(tex "\\begin{split}\nu_{i,j}^{n+1} = & u_{i,j}^n - \\frac{\\Delta t}{\\Delta x} u_{i,j}^n (u_{i,j}^n - u_{i-1,j}^n)  - \\frac{\\Delta t}{\\Delta y} v_{i,j}^n (u_{i,j}^n - u_{i,j-1}^n) \\\\\n&+ \\frac{\\nu \\Delta t}{\\Delta x^2}(u_{i+1,j}^n-2u_{i,j}^n+u_{i-1,j}^n) + \\frac{\\nu \\Delta t}{\\Delta y^2} (u_{i,j+1}^n - 2u_{i,j}^n + u_{i,j-1}^n)\n\\end{split}")
;;
(tex "\\begin{split}\nv_{i,j}^{n+1} = & v_{i,j}^n - \\frac{\\Delta t}{\\Delta x} u_{i,j}^n (v_{i,j}^n - v_{i-1,j}^n) - \\frac{\\Delta t}{\\Delta y} v_{i,j}^n (v_{i,j}^n - v_{i,j-1}^n) \\\\\n&+ \\frac{\\nu \\Delta t}{\\Delta x^2}(v_{i+1,j}^n-2v_{i,j}^n+v_{i-1,j}^n) + \\frac{\\nu \\Delta t}{\\Delta y^2} (v_{i,j+1}^n - 2v_{i,j}^n + v_{i,j-1}^n)\n\\end{split}")
;;
;; ### Initial Conditions
;;
(def nx 41)
(def ny 41)
(def nt 120)
(def c 1)
(def nu 0.01)
(def sigma 0.0009)
(def dx (/ 2 (- nx 1)))
(def dy (/ 2 (- ny 1)))
(def init-params
  {:nx    nx
   :ny    ny
   :nt    nt
   :c     c
   :nu    nu
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
;; Create the initial u and v arrays w/ initial conditions:
;; u(.5<=x<=1 && .5<=y<=1) is 2, else 1
;; v(.5<=x<=1 && .5<=y<=1) is 2, else 1
;;
(def array-u (two-d/create-init-u init-params spatial-arr))
(def array-v (two-d/create-init-u init-params spatial-arr))
;;
;; initial plotting w/ $x$, $y$, $u$ goes:
;;
(two-d/sim->plotly-plot-it! spatial-arr array-u)
;;
(def sim-result (two-d/simulate {:array-u array-u
                                 :array-v array-v} (assoc init-params :mode :burgers)))
;;
^:kindly/hide-code
(let [{:keys [array-u array-v]} sim-result
      plottable-data-1 (two-d/arr->plotly-plottable-data spatial-arr array-u)
      plottable-data-2 (two-d/arr->plotly-plottable-data spatial-arr array-v)]
  (kind/plotly
    {:data   [(merge two-d/plotly-opts plottable-data-1)
              (merge two-d/plotly-opts plottable-data-2)]
     :layout {:scene {:zaxis {:range [0.8 2.2]}}}}))
;;
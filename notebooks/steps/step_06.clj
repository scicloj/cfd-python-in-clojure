^:kindly/hide-code
(ns steps.step-06
  (:require
   [cfd.two-d :as two-d]
   [scicloj.kindly.v4.api :as kindly]
   [scicloj.kindly.v4.kind :as kind]
   [tablecloth.api :as tc]
   [utils.notebook :refer [md tex]]))

;; ## Step 6: 2-D Convection
;;
;; Now we solve 2D Convection, represented by the pair of coupled partial differential equations below:
;;
(tex "\\frac{\\partial u}{\\partial t} + u \\frac{\\partial u}{\\partial x} + v \\frac{\\partial u}{\\partial y} = 0")
;;
(tex "\\frac{\\partial v}{\\partial t} + u \\frac{\\partial v}{\\partial x} + v \\frac{\\partial v}{\\partial y} = 0")
;;
;; Descretizing these equations using hte methods we've applied previously yields:
;;
(tex "\\frac{u_{i,j}^{n+1}-u_{i,j}^n}{\\Delta t} + u_{i,j}^n \\frac{u_{i,j}^n-u_{i-1,j}^n}{\\Delta x} + v_{i,j}^n \\frac{u_{i,j}^n-u_{i,j-1}^n}{\\Delta y} = 0")
;;
(tex "\\frac{v_{i,j}^{n+1}-v_{i,j}^n}{\\Delta t} + u_{i,j}^n \\frac{v_{i,j}^n-v_{i-1,j}^n}{\\Delta x} + v_{i,j}^n \\frac{v_{i,j}^n-v_{i,j-1}^n}{\\Delta y} = 0")
;;
;; Rearranging both equations, we solve for $u_{i,j}^{n+1}$ and $v_{i,j}^{n+1}$, respectively.
;; Note that these equations are also coupled.
;;
(tex "u_{i,j}^{n+1} = u_{i,j}^n - u_{i,j} \\frac{\\Delta t}{\\Delta x} (u_{i,j}^n-u_{i-1,j}^n) - v_{i,j}^n \\frac{\\Delta t}{\\Delta y} (u_{i,j}^n-u_{i,j-1}^n)")
;;
(tex "v_{i,j}^{n+1} = v_{i,j}^n - u_{i,j} \\frac{\\Delta t}{\\Delta x} (v_{i,j}^n-v_{i-1,j}^n) - v_{i,j}^n \\frac{\\Delta t}{\\Delta y} (v_{i,j}^n-v_{i,j-1}^n)")
;;
;; ### Initial Conditions
;;
;; The initial conditions are the same that we used for 1D convection, applied in both the $x$ and $y$ directions.
;;
(tex "u,\\ v\\ = \\begin{cases}\\begin{matrix}\n2 & \\text{for } x,y \\in (0.5, 1)\\times(0.5,1) \\cr\n1 & \\text{everywhere else}\n\\end{matrix}\\end{cases}")
;;
;; ### Boundary Conditions
;;
;; The boundary conditions hold $u$ and $v$ equal to 1 along the boundaries of the grid.
;;
(tex "u = 1,\\ v = 1 \\text{ for } \\begin{cases} \\begin{matrix}x=0,2\\cr y=0,2 \\end{matrix}\\end{cases}")
;;
;; number of x grid points
(def nx 101)
;; number of y grid points
(def ny 101)
;; number of time steps
(def nt 80)
;; CFL number
(def sigma 0.2)
;; Define initial parameters
(def init-params
  {:nx    nx
   :ny    ny
   :nt    nt
   :c     1
   :dx    (/ 2 (- nx 1))
   :dy    (/ 2 (- ny 1))
   :sigma sigma
   :dt    (* sigma (/ 2 (- nx 1)))})
;;
;; Create the spatial grid
(def grid-start 0) (def grid-end 2)
(def spatial-arr (two-d/create-array-2d
                   (assoc init-params
                     :x-start grid-start :x-end grid-end
                     :y-start grid-start :y-end grid-end)))
;;
;; Create the initial u and v arrays
(def array-u (two-d/create-init-u init-params spatial-arr))
(def array-v (two-d/create-init-u init-params spatial-arr))
;;
^:kindly/hide-code
(defn arr->plotly-plottable-data
  [vel-array [array-x array-y]]
  {:x (apply concat (repeat (alength array-y) array-x))
   :y (apply concat (map #(repeat (alength array-x) %) array-y))
   :z (apply concat vel-array)})
;;
^:kindly/hide-code
(def plotly-plottable-data (arr->plotly-plottable-data array-u spatial-arr))
;;
^:kindly/hide-code
(-> plotly-plottable-data
    tc/dataset
    (kind/dataset {:dataset/print-range 6}))
;;
;; initial plotting w/ $x$, $y$, $u$ goes:
;;
^:kindly/hide-code
(def plotly-opts {:type    :mesh3d
                  :opacity 0.20
                  :color   "blue"
                  :marker  {:colorscale :Viridis}})
;;
^:kindly/hide-code
(kind/plotly {:data [(merge plotly-plottable-data plotly-opts)]})

;; ### Iterating in 2-D w/ nonlinear convection equation
;;
^:kindly/hide-code
(def simulated-result (two-d/simulate
                        {:array-u array-u
                         :array-v array-v}
                        (assoc init-params :mode :nonlinear-convection)))


;; u
(let [plottable-data (arr->plotly-plottable-data (:array-u simulated-result) spatial-arr)]
  (kind/plotly
    {:data [(merge plottable-data plotly-opts)]}))
;; v
(let [plottable-data (arr->plotly-plottable-data (:array-v simulated-result) spatial-arr)]
  (kind/plotly
    {:data [(merge plottable-data (assoc plotly-opts :color "red"))]}))
;;
;;
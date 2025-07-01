^:kindly/hide-code
(ns steps.step-05
  (:require
   [cfd.two-d :as two-d]
   [scicloj.kindly.v4.kind :as kind]
   [tablecloth.api :as tc]
   [utils.notebook :refer [md tex]]))

;; Expanding from 1D to 2D.
;; Following exercises extends firstly to 2D. The expansion simply
;; requires to apply the definition: a partial derivative with respect
;; to $x$ is the variation in the $x$ direction at constant $y$.

;; In 2D space, a rectangular(uniform) grid is defined by the points with coordinates:

(tex "x_i = x_0 + i \\Delta x")
(tex "y_i = y_0 + i \\Delta y")

;; Then also define $u_{i,j} = u(x_i, y_j) and apply the finite-difference formulas
;; on either variable $x, y$ _acting separately_ on the $i$ and $j$ indices.
;; All derivatives are based on the 2D Taylor expansion of a mesh point value around $u_{ij}$
;;
;; Hence, for a first-order partial derivative in the x-direction, a finite-difference formula is:
;;
(tex "\\frac{\\partial u}{\\partial x}\\biggr\\rvert_{i,j} = \\frac{u_{i+1,j}-u_{i,j}}{\\Delta x}+\\mathcal{O}(\\Delta x)")
;;
;; and similarly in the $y$ direction. Thus, we can write backward-difference,
;; forward-difference or central difference formulas for Step 5 to 12.
;;
;; ## Step 5: 2-D Linear Convection
;;
;; The PDE governing 2-D Linear Convection is written as
(tex "\\frac{\\partial u}{\\partial t}+c\\frac{\\partial u}{\\partial x} + c\\frac{\\partial u}{\\partial y} = 0")
;; This is the same form in 1-D, then added one more dimension to account for
;; as we step forward in time.

;; We will use:
;; - **a forward difference discretization** for the timestep
;; - **a backward difference discretization** for two spatial steps
;;
;; With 1-D implementations, we used $i$ subscripts to denote movement in space
;; (e.g. $u_i^n - u_{i-1}^n$). Now that we have two dimensions to account for,
;; we need to add a second subscript, $j$, to account for all the information in
;; the regime.
;;
;; Here, we'll again use $i$ as the index for our $x$ values, and we'll add
;; the $j$ subscript to track our $y$ values.
;;
;; With that in mind, our discretization of the PD should be relatively straightforward.
;;
(tex "\\frac{u_{i,j}^{n+1}-u_{i,j}^n}{\\Delta t} + c\\frac{u_{i, j}^n-u_{i-1,j}^n}{\\Delta x} + c\\frac{u_{i,j}^n-u_{i,j-1}^n}{\\Delta y}=0")
;;
;; As before, solve for the only unknown:
;;
(tex "u_{i,j}^{n+1} = u_{i,j}^n-c \\frac{\\Delta t}{\\Delta x}(u_{i,j}^n-u_{i-1,j}^n)-c \\frac{\\Delta t}{\\Delta y}(u_{i,j}^n-u_{i,j-1}^n)")
;;
;; We will solve this equation with the following initial conditions:
;;
(tex "u(x,y) = \\begin{cases}\n\\begin{matrix}\n2\\ \\text{for} & 0.5 \\leq x, y \\leq 1 \\cr\n1\\ \\text{for} & \\text{everywhere else}\\end{matrix}\\end{cases}")
;;
;; and the boundary conditions:
;;
(tex "u = 1\\ \\text{for } \\begin{cases}\n\\begin{matrix}\nx =  0,\\ 2 \\cr\ny =  0,\\ 2 \\end{matrix}\\end{cases}")
;;
^:kindly/hide-code (def nx 81)
^:kindly/hide-code (def ny 81)
^:kindly/hide-code (def sigma 0.2)
^:kindly/hide-code (def dx (/ 2 (- nx 1)))
^:kindly/hide-code
(def init-params
  {:nx    nx
   :ny    ny
   :nt    100
   :c     1
   :dx    dx
   :dy    (/ 2 (- ny 1))
   :sigma sigma
   :dt    (* sigma dx)})

^:kindly/hide-code
(def spatial-arr (two-d/create-2d-spacial-array
                   {:nx nx :x-start 0 :x-end 2
                    :ny ny :y-start 0 :y-end 2}))

^:kindly/hide-code
(def array-u (two-d/create-init-u (merge init-params spatial-arr)))

;; We plot here using plotly, then using `:mesh3d` as the type of the plot.
;; And here's [a reference doc](https://scicloj.github.io/kindly-noted/kinds.html#plotly) from kindly notebook.
;; The plotting data formats goes like:

^:kindly/hide-code
(defn arr->plotly-plottable-data
  [{:keys [array-x array-y array-u]}]
  {:x (apply concat (repeat (alength array-y) array-x))
   :y (apply concat (map #(repeat (alength array-x) %) array-y))
   :z (apply concat array-u)})

^:kindly/hide-code
(def plotly-plottable-data (arr->plotly-plottable-data (assoc spatial-arr :array-u array-u)))

^:kindly/hide-code
(-> plotly-plottable-data
    tc/dataset
    (kind/dataset {:dataset/print-range 6}))

;;
;; initial plotting goes:
;;
^:kindly/hide-code
(def plotly-opts {:type    :mesh3d
                  :opacity 0.20
                  :color   "lightpink"
                  :marker  {:colorscale :Viridis}})
^:kindly/hide-code
(kind/plotly
  {:data [(merge plotly-plottable-data plotly-opts)]})

;;
;; **note**: for now, we skip 3d plotting notes from PythonCFD(further **_TODO_**)
;;

;; ### Iterating in 2-D w/ linear convection equation
;;
;;
^:kindly/hide-code
(let [convection-u   (two-d/simulate array-u init-params)
      plottable-data (arr->plotly-plottable-data (assoc spatial-arr :array-u convection-u))]
  (kind/plotly
    {:data [(merge plottable-data plotly-opts)]}))

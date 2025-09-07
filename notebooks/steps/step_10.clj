^:kindly/hide-code
(ns steps.step-10
  (:require
   [cfd.two-d :as two-d]
   [fastmath.core :as fm :refer [pow]]
   [scicloj.kindly.v4.api :as kindly]
   [scicloj.kindly.v4.kind :as kind]
   [utils.notebook :refer [md tex]]))

;; For a moment, recall the Navier-Stokes equations for an incompressible fluid,
;; where $\vec{v}$ represents the velocity field:
(tex "\\begin{eqnarray*}\n\\nabla \\cdot\\vec{v} &=& 0 \\\\\n\\frac{\\partial \\vec{v}}{\\partial t}+(\\vec{v}\\cdot\\nabla)\\vec{v} &=& -\\frac{1}{\\rho}\\nabla p + \\nu \\nabla^2\\vec{v}\n\\end{eqnarray*}")
;; The first equation represents mass conservation at constant density.
;; The second equation is the conservation of momentum. But a problem appears:
;; the continuity equation for incompressible flow does not have a dominant variable, and
;; there is no obvious way to couple the velocity and the pressure.
;; In the case of compressible flow, in contrast, mass continuity would provide an evolution
;; for the density $\rho$, which is coupled with an equation of state relating $\rho$ and $p$.
;;
;; In incompressible flow, the continuity equation $\nabla \cdot\vec{v}=0$ provides a
;; _kinematic constraint_ that requires the pressure field to evolve so that the rate of expansion
;; $\nabla \cdot\vec{v}$ should vanish everywhere. A way out of this difficulty is to _construct_
;; a pressure field that guarantees continuity is satisfied; such a relation can be obtained by
;; taking the divergence of the momentum equation. In that process, a Poisson equation
;; for the pressure shows up.
;;
;; ## Step 10: 2D Poisson Equation
;;
;; Poisson's equation is obtained from adding a source term to the right-hand-side of Laplace's equation:
(tex "\\frac{\\partial ^2 p}{\\partial x^2} + \\frac{\\partial ^2 p}{\\partial y^2} = b")
;;
;; So, unlike the Laplace equation, there is some finite value inside the field that affects
;; the solution. Poisson's equation acts to "relax" the initial source in the field.
;;
;; In discretized form, this looks almost the same as [Step 9](/steps.step_09.html), except for the source term:
(tex "\\frac{p_{i+1,j}^{n}-2p_{i,j}^{n}+p_{i-1,j}^{n}}{\\Delta x^2}+\\frac{p_{i,j+1}^{n}-2 p_{i,j}^{n}+p_{i,j-1}^{n}}{\\Delta y^2}=b_{i,j}^{n}")
;;
;; As before, we rearrange this so that we obtain an equation for $p$ at point $i, j$.
;; Thus, we obtain:
(tex "p_{i,j}^{n}=\\frac{(p_{i+1,j}^{n}+p_{i-1,j}^{n})\\Delta y^2+(p_{i,j+1}^{n}+p_{i,j-1}^{n})\\Delta x^2-b_{i,j}^{n}\\Delta x^2\\Delta y^2}{2(\\Delta x^2+\\Delta y^2)}")
;;
;; We will solve this equation by assuming an initial state of $p = 0$ everywhere, and
;; applying boundary conditions as follows:
;;
;; $p = 0$ at $x = 0, 2$ and $y = 0, 1$
;;
;; and the source term consists of two initial spikes inside the domain, as follows:
;;
;; $b_{i,j}=100$ at $i=\frac{1}{4}nx, j=\frac{1}{4}ny$
;;
;; $b_{i,j}=-100$ at $i=\frac{3}{4}nx, j=\frac{3}{4}ny$
;;
;; $b_{i,j}=0$ everywhere else.
;;
;; The iterations will advance in pseudo-time to relax the initial spikes. The relaxation under
;; Poisson's equation gets slower and slower as they pregress. _Why?_
;;
(def nx 50)
(def ny 50)
(def nt 100)
(def nt 100)
(def dx (double (/ 2 (dec nx))))
(def dy (double (/ 1 (dec ny))))
(def spatial-init-param
  {:nx nx :x-start 0 :x-end 2 :dx dx :dx-square (* dx dx)
   :ny ny :y-start 0 :y-end 1 :dy dy :dy-square (* dy dy)
   :nt nt})
;;
;; Have the discretized 2D spatial array ready:
(def spatial-array (two-d/create-array-2d spatial-init-param))

(def array-p (two-d/create-init-u
               (assoc spatial-init-param
                      :d-type Double/TYPE
                      :condition-fn (constantly 0.0))
               spatial-array))

(def array-b (aclone array-p))
(aset array-b (int (* ny 0.25)) (int (* nx 0.25)) 100.0)
(aset array-b (int (* ny 0.25 3)) (int (* nx 0.25 3)) -100.0)

;; With the initial parameters setup above, we are ready to advance
;; the initial guess in pseudo-time. How is the code below difference
;; from the function used in [Step 9](/steps.step_09.html) to solve
;; Laplace's equation?

(defn poisson-2d [{:keys [spatial-array
                          array-p
                          array-b
                          nx ny
                          dx dy
                          dx-square dy-square]
                   :as   params}]
  (let [pd (two-d/clone-2d-array array-p)]
   (dotimes [y-idx (- ny 2)]
     (let [prev-y-idx y-idx
           y-idx      (inc y-idx)
           next-y-idx (inc y-idx)

           pd-prev-x  (aget pd prev-y-idx)
           pd-x       (aget pd y-idx)
           pd-next-x  (aget pd next-y-idx)]
      (dotimes [x-idx (- nx 2)]
        (let [prev-x-idx x-idx
              x-idx      (inc x-idx)
              next-x-idx (inc x-idx)

              p-j-i+1    (aget pd-x next-x-idx)
              p-j-i-1    (aget pd-x prev-x-idx)
              p-j+1-i    (aget pd-next-x x-idx)
              p-j-1-i    (aget pd-prev-x x-idx)

              b-j-i      (aget array-b y-idx x-idx)]
          (aset array-p y-idx x-idx
            (/ (+ (* dy-square (+ p-j-i+1 p-j-i-1))
                  (* dx-square (+ p-j+1-i p-j-1-i))
                  (* -1.0 b-j-i dx-square dy-square))
               (* 2 (+ dx-square dy-square)))))))))

  ;; Boundary conditions
  (dotimes [y-idx ny]
    (aset array-p y-idx 0 0.0)
    (aset array-p y-idx (dec nx) 0.0))
  (dotimes [x-idx nx]
    (aset array-p 0 x-idx 0.0)
    (aset array-p (dec ny) x-idx 0.0)))

(defn simulate [{:keys [array-p] :as params}]
  (loop [n 0]
    (if (= n nt)
      array-p
      (do
        (poisson-2d params)
        (recur (inc n))))))

(def init-params
  (assoc spatial-init-param
         :array-p array-p
         :array-b array-b))

(def plotly-surface-plot-base-opts
  (let [[array-x array-y] spatial-array]
    {:layout {:scene {:xaxis {:range [0.0 2.1]}
                      :yaxis {:range [0.0 1.1]}
                      :zaxis {:range [-0.2 0.2]}}}
     :data   [{:x       array-x
               :y       array-y
               :z       array-p
               :type    :surface
               :opacity 0.20
               :color   "size"
               :marker  {:colorscale :Viridis}}]}))

;; Run simulation and plot the result:
(let [simulated-data (simulate init-params)]
  (kind/plotly
    (update-in plotly-surface-plot-base-opts
      [:data 0] #(assoc % :z simulated-data))))

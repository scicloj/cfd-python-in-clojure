^:kindly/hide-code
(ns steps.step-11
  (:require
   [cfd.two-d :as two-d]
   [fastmath.core :as fm :refer [pow]]
   [scicloj.kindly.v4.api :as kindly]
   [scicloj.kindly.v4.kind :as kind]
   [utils.notebook :refer [md tex]]))

;; The final two steps in this interactive module porting CFD Python into Clojure will both solve
;; both the Navier-Stokes equations in two dimensions, but with different boundary conditions.
;;
;; The momentum equation in vector form for a velocity field $\vec{v}$is:
(tex "\\frac{\\partial \\vec{v}}{\\partial t}+(\\vec{v}\\cdot\\nabla)\\vec{v}=-\\frac{1}{\\rho}\\nabla p + \\nu \\nabla^2\\vec{v}")
;;
;; This represents three scalar equations, one for each velocity component $(u, v, w)$. But we will
;; solve it in two dimensions, so there will be two scalar equations.
;;
;; Remember the continuity equation? This is where the [Poisson equation](/steps.step_10.html) for pressure comes in!
;;
;; ## Step 11: Cavity Flow with Navier-Stokes
;;
;; Here is the system of differential equations: two equations for the velocity components $u$,
;; $v$ and one equation for pressure:
(tex "\\frac{\\partial u}{\\partial t}+u\\frac{\\partial u}{\\partial x}+v\\frac{\\partial u}{\\partial y} = -\\frac{1}{\\rho}\\frac{\\partial p}{\\partial x}+\\nu \\left(\\frac{\\partial^2 u}{\\partial x^2}+\\frac{\\partial^2 u}{\\partial y^2} \\right)")
;;
(tex "\\frac{\\partial v}{\\partial t}+u\\frac{\\partial v}{\\partial x}+v\\frac{\\partial v}{\\partial y} = -\\frac{1}{\\rho}\\frac{\\partial p}{\\partial y}+\\nu\\left(\\frac{\\partial^2 v}{\\partial x^2}+\\frac{\\partial^2 v}{\\partial y^2}\\right)")
;;
(tex "\\frac{\\partial^2 p}{\\partial x^2}+\\frac{\\partial^2 p}{\\partial y^2} = -\\rho\\left(\\frac{\\partial u}{\\partial x}\\frac{\\partial u}{\\partial x}+2\\frac{\\partial u}{\\partial y}\\frac{\\partial v}{\\partial x}+\\frac{\\partial v}{\\partial y}\\frac{\\partial v}{\\partial y} \\right)")
;;
;; From the previous steps, we already know how to discretize all these terms. Only the last equation
;; is a little unfamiliar. But with a little patience, it will not be hard!
;;
;; ### Discretized equations
;;
;; First, let's discretize the $u$-momentum equation, as follows:
(tex "\\begin{split}\n& \\frac{u_{i,j}^{n+1}-u_{i,j}^{n}}{\\Delta t}+u_{i,j}^{n}\\frac{u_{i,j}^{n}-u_{i-1,j}^{n}}{\\Delta x}+v_{i,j}^{n}\\frac{u_{i,j}^{n}-u_{i,j-1}^{n}}{\\Delta y} = \\\\ \n& \\qquad -\\frac{1}{\\rho}\\frac{p_{i+1,j}^{n}-p_{i-1,j}^{n}}{2\\Delta x}+\\nu\\left(\\frac{u_{i+1,j}^{n}-2u_{i,j}^{n}+u_{i-1,j}^{n}}{\\Delta x^2}+\\frac{u_{i,j+1}^{n}-2u_{i,j}^{n}+u_{i,j-1}^{n}}{\\Delta y^2}\\right)\n\\end{split}")
;;
;; Similarly for the $v$-momentum equation:
(tex "\\begin{split}\n&\\frac{v_{i,j}^{n+1}-v_{i,j}^{n}}{\\Delta t}+u_{i,j}^{n}\\frac{v_{i,j}^{n}-v_{i-1,j}^{n}}{\\Delta x}+v_{i,j}^{n}\\frac{v_{i,j}^{n}-v_{i,j-1}^{n}}{\\Delta y} = \\\\\n& \\qquad -\\frac{1}{\\rho}\\frac{p_{i,j+1}^{n}-p_{i,j-1}^{n}}{2\\Delta y}\n+\\nu\\left(\\frac{v_{i+1,j}^{n}-2v_{i,j}^{n}+v_{i-1,j}^{n}}{\\Delta x^2}+\\frac{v_{i,j+1}^{n}-2v_{i,j}^{n}+v_{i,j-1}^{n}}{\\Delta y^2}\\right)\n\\end{split}")
;;
;; Finally, the discretized pressure-Poisson equation can be written thus:
(tex "\\begin{split}\n& \\frac{p_{i+1,j}^{n}-2p_{i,j}^{n}+p_{i-1,j}^{n}}{\\Delta x^2}+\\frac{p_{i,j+1}^{n}-2p_{i,j}^{n}+p_{i,j-1}^{n}}{\\Delta y^2} = \\\\\n& \\qquad \\rho \\left[ \\frac{1}{\\Delta t}\\left(\\frac{u_{i+1,j}-u_{i-1,j}}{2\\Delta x}+\\frac{v_{i,j+1}-v_{i,j-1}}{2\\Delta y}\\right) -\\frac{u_{i+1,j}-u_{i-1,j}}{2\\Delta x}\\frac{u_{i+1,j}-u_{i-1,j}}{2\\Delta x} - 2\\frac{u_{i,j+1}-u_{i,j-1}}{2\\Delta y}\\frac{v_{i+1,j}-v_{i-1,j}}{2\\Delta x} - \\frac{v_{i,j+1}-v_{i,j-1}}{2\\Delta y}\\frac{v_{i,j+1}-v_{i,j-1}}{2\\Delta y}\\right]\n\\end{split}")
;;
;; You should write these equations down on your own notes, by hand, following each term mentally as you write it.
;;
;; As before, let's rearrange the equations in the way that the iterations need to proceed in the code.
;; First, the momentum equations for the velocity at the next time step.
;;
;; The momentum equation in the $u$ direction:
(tex "\\begin{split}\nu_{i,j}^{n+1} = u_{i,j}^{n} & - u_{i,j}^{n} \\frac{\\Delta t}{\\Delta x} \\left(u_{i,j}^{n}-u_{i-1,j}^{n}\\right) - v_{i,j}^{n} \\frac{\\Delta t}{\\Delta y} \\left(u_{i,j}^{n}-u_{i,j-1}^{n}\\right) \\\\\n& - \\frac{\\Delta t}{\\rho 2\\Delta x} \\left(p_{i+1,j}^{n}-p_{i-1,j}^{n}\\right) \\\\\n& + \\nu \\left(\\frac{\\Delta t}{\\Delta x^2} \\left(u_{i+1,j}^{n}-2u_{i,j}^{n}+u_{i-1,j}^{n}\\right) + \\frac{\\Delta t}{\\Delta y^2} \\left(u_{i,j+1}^{n}-2u_{i,j}^{n}+u_{i,j-1}^{n}\\right)\\right)\n\\end{split}")
;;;
;; The momentum equation in the $v$ direction:
(tex "\\begin{split}\nv_{i,j}^{n+1} = v_{i,j}^{n} & - u_{i,j}^{n} \\frac{\\Delta t}{\\Delta x} \\left(v_{i,j}^{n}-v_{i-1,j}^{n}\\right) - v_{i,j}^{n} \\frac{\\Delta t}{\\Delta y} \\left(v_{i,j}^{n}-v_{i,j-1}^{n})\\right) \\\\\n& - \\frac{\\Delta t}{\\rho 2\\Delta y} \\left(p_{i,j+1}^{n}-p_{i,j-1}^{n}\\right) \\\\\n& + \\nu \\left(\\frac{\\Delta t}{\\Delta x^2} \\left(v_{i+1,j}^{n}-2v_{i,j}^{n}+v_{i-1,j}^{n}\\right) + \\frac{\\Delta t}{\\Delta y^2} \\left(v_{i,j+1}^{n}-2v_{i,j}^{n}+v_{i,j-1}^{n}\\right)\\right)\n\\end{split}")
;;
;; Almost there! Now, we rearrange the pressure-Poisson equation:
(tex "\\begin{split}\np_{i,j}^{n} = & \\frac{\\left(p_{i+1,j}^{n}+p_{i-1,j}^{n}\\right) \\Delta y^2 + \\left(p_{i,j+1}^{n}+p_{i,j-1}^{n}\\right) \\Delta x^2}{2\\left(\\Delta x^2+\\Delta y^2\\right)} \\\\\n& -\\frac{\\rho\\Delta x^2\\Delta y^2}{2\\left(\\Delta x^2+\\Delta y^2\\right)} \\\\\n& \\times \\left[\\frac{1}{\\Delta t}\\left(\\frac{u_{i+1,j}-u_{i-1,j}}{2\\Delta x}+\\frac{v_{i,j+1}-v_{i,j-1}}{2\\Delta y}\\right)-\\frac{u_{i+1,j}-u_{i-1,j}}{2\\Delta x}\\frac{u_{i+1,j}-u_{i-1,j}}{2\\Delta x} -2\\frac{u_{i,j+1}-u_{i,j-1}}{2\\Delta y}\\frac{v_{i+1,j}-v_{i-1,j}}{2\\Delta x}-\\frac{v_{i,j+1}-v_{i,j-1}}{2\\Delta y}\\frac{v_{i,j+1}-v_{i,j-1}}{2\\Delta y}\\right]\n\\end{split}")
;;
;; The initial condition is $u, v, p = 0$ everywhere, and the boundary conditions are:
;;
;; $u = 1$ at $y = 2$ (the "lid") <br>
;; $u, v = 0$ on the other boundaries <br>
;; $\frac{\partial p}{\partial y} = 0$ at $y = 0$ <br>
;; $p = 0$ at $y = 2$ <br>
;; $\frac{\partial p}{\partial x} = 0$ at $x = 0, 2$
;;
;; ### Implementing Cavity Flow
;;
;; Spatial setup:
(def nx 41)
(def ny 41)
(def nt 500)
(def nit 50)
(def c 1)
(def x-start 0)
(def x-end 2)
(def y-start 0)
(def y-end 2)
(def spatial-init-params
  {:nx nx :x-start x-start :x-end x-end :dx (double (/ (- x-end x-start) (dec nx)))
   :ny ny :y-start y-start :y-end y-end :dy (double (/ (- y-end y-start) (dec ny)))})

(def spatial-array (two-d/create-array-2d spatial-init-params))
;;
;; The rest setup:
(def rho 1.0)
(def nu 0.1)
(def dt 0.001)

(def array-u (two-d/create-init-zeros-u (assoc spatial-init-params :d-type Double/TYPE) spatial-array))
(def array-v (two-d/clone-2d-array array-u))
(def array-p (two-d/clone-2d-array array-u))
(def array-b (two-d/clone-2d-array array-u))

(def init-params
  (assoc spatial-init-params
         :rho rho
         :nu nu
         :dt dt
         :nit nit
         :nt nt
         :spatial-array spatial-array
         :array-u array-u
         :array-v array-v
         :array-p array-p
         :array-b array-b))
;;
;; The pressure Poisson equation that's written above can be hard to write out without typos.
;; The function `build-up-b` below represents the contents of the square brackets, so that the entirety of
;; the PPE is slightly more manageable.
;;
(defn build-up-b [{:keys [array-b array-u array-v rho dt dx dy nx ny] :as _params}]
  (dotimes [y-idx (- ny 2)]
    (dotimes [x-idx (- nx 2)]
      (let [y-idx   (inc y-idx)
            x-idx   (inc x-idx)
            u-i-j-1 (aget array-u (dec y-idx) x-idx)
            u-i-j+1 (aget array-u (inc y-idx) x-idx)
            u-j-i-1 (aget array-u y-idx (dec x-idx))
            u-j-i+1 (aget array-u y-idx (inc x-idx))
            v-i-j-1 (aget array-v (dec y-idx) x-idx)
            v-i-j+1 (aget array-v (inc y-idx) x-idx)
            v-j-i-1 (aget array-v y-idx (dec x-idx))
            v-j-i+1 (aget array-v y-idx (inc x-idx))]
        (aset array-b y-idx x-idx
          (* rho
             (- (/ (+ (/ (- u-j-i+1 u-j-i-1)
                         (* 2.0 dx))
                      (/ (- v-i-j+1 v-i-j-1)
                         (* 2.0 dy)))
                   dt)
                (/ (pow (- u-j-i+1 u-j-i-1) 2)
                   (pow (* 2.0 dx) 2))
                (/ (pow (- v-i-j+1 v-i-j-1) 2)
                   (pow (* 2.0 dy) 2))
                (/ (* (- u-i-j+1 u-i-j-1) (- v-j-i+1 v-j-i-1))
                   (* 2.0 dx dy))))))))
  array-b)

;; The function `pressure-poisson` is also defined to help segregate the different rounds of calculations.
;; Note the presence of the pseudo-time variable `nit`.
;; This sub-iteration in the Poisson calculation helps
;; ensure a divergence-free field.
;;
(defn pressure-poisson [{:keys [array-p array-b nx ny dx dy nit] :as _params}]
  (loop [t-idx 0]
    (if (= t-idx nit)
      array-p
      (let [pn (two-d/clone-2d-array array-p)]
        (dotimes [y-idx (- ny 2)]
          (dotimes [x-idx (- nx 2)]
            (let [y-idx   (inc y-idx)
                  x-idx   (inc x-idx)
                  p-i-j-1 (aget pn (dec y-idx) x-idx)
                  p-i-j+1 (aget pn (inc y-idx) x-idx)
                  p-j-i-1 (aget pn y-idx (dec x-idx))
                  p-j-i+1 (aget pn y-idx (inc x-idx))]
              (aset array-p y-idx x-idx
                (- (/ (+ (* (+ p-j-i+1 p-j-i-1) (pow dy 2))
                         (* (+ p-i-j+1 p-i-j-1) (pow dx 2)))
                      (* 2.0 (+ (pow dx 2) (pow dy 2))))
                   (* (/ (* (pow dx 2) (pow dy 2))
                         (* 2.0 (+ (pow dx 2) (pow dy 2))))
                      (aget array-b y-idx x-idx)))))))
        ;; boundary conditions
        ;; dp/dx = 0 at x = 2 & dp/dx = 0 at x = 0
        (dotimes [y-idx ny]
          (aset array-p y-idx (dec nx) (aget array-p y-idx (- nx 2)))
          (aset array-p y-idx 0 (aget array-p y-idx 1)))
        ;; dp/dy = 0 at y = 0 & p = 0 at y = 2
        (dotimes [x-idx nx]
          (aset array-p 0 x-idx (aget array-p 1 x-idx))
          (aset array-p (dec ny) x-idx 0.0))
        (recur (inc t-idx))))))

;; Finally, the rest of the cavity flow equations are wrapped inside the function `cavity-flow`,
;; allowing us to easily plot the results of the cavity flow solver for different lengths of time.
;;
(defn cavity-flow [{:keys [nt array-u array-v array-p dt nx ny dx dy rho nu] :as params}]
  (loop [n 0]
    (if (= n nt)
      params
      (let [un (two-d/clone-2d-array array-u)
            vn (two-d/clone-2d-array array-v)]
        (build-up-b params)
        (pressure-poisson params)
        (dotimes [y-idx (- ny 2)]
          (dotimes [x-idx (- nx 2)]
            (let [y-idx   (inc y-idx)
                  x-idx   (inc x-idx)
                  u-i-j   (aget un y-idx x-idx)
                  u-i-j-1 (aget un (dec y-idx) x-idx)
                  u-i-j+1 (aget un (inc y-idx) x-idx)
                  u-j-i-1 (aget un y-idx (dec x-idx))
                  u-j-i+1 (aget un y-idx (inc x-idx))
                  v-i-j   (aget vn y-idx x-idx)
                  v-i-j-1 (aget vn (dec y-idx) x-idx)
                  v-i-j+1 (aget vn (inc y-idx) x-idx)
                  v-j-i-1 (aget vn y-idx (dec x-idx))
                  v-j-i+1 (aget vn y-idx (inc x-idx))
                  p-j-i+1 (aget array-p y-idx (inc x-idx))
                  p-j-i-1 (aget array-p y-idx (dec x-idx))]
              (aset array-u y-idx x-idx
                (- u-i-j
                   (* u-i-j dt (/ 1.0 dx) (- u-i-j u-j-i-1))
                   (* v-i-j dt (/ 1.0 dy) (- u-i-j u-j-i-1))
                   (* dt (/ 1.0 (* 2.0 rho dx)) (- p-j-i+1 p-j-i-1))
                   (* -1.0
                      nu
                      (+ (* dt
                            (/ 1.0 (pow dx 2))
                            (+ u-j-i+1
                               (* -2.0 u-i-j)
                               u-j-i-1))
                         (* dt
                            (/ 1.0 (pow dy 2))
                            (+ u-i-j+1
                               (* -2.0 u-i-j)
                               u-i-j-1))))))

              (aset array-v y-idx x-idx
                (- v-i-j
                   (* v-i-j dt (/ 1.0 dx) (- v-i-j v-j-i-1))
                   (* u-i-j dt (/ 1.0 dy) (- v-i-j v-j-i-1))
                   (* dt (/ 1.0 (* 2.0 rho dy)) (- p-j-i+1 p-j-i-1))
                   (* -1.0
                      nu
                      (+ (* dt
                            (/ 1.0 (pow dx 2))
                            (+ v-j-i+1
                               (* -2.0 v-i-j)
                               v-j-i-1))
                         (* dt
                            (/ 1.0 (pow dy 2))
                            (+ v-i-j+1
                               (* -2.0 v-i-j)
                               v-i-j-1)))))))))

        ;; boundary conditions
        (dotimes [y-idx ny]
          (aset array-u y-idx 0 0.0)
          (aset array-u y-idx (dec nx) 1.0)                 ;; set velocity on cavity lid equal to 1
          (aset array-v y-idx 0 0.0)
          (aset array-v y-idx (dec nx) 0.0))
        (dotimes [x-idx nx]
          (aset array-u 0 x-idx 0.0)
          (aset array-u (dec ny) x-idx 0.0)
          (aset array-u 0 x-idx 0.0)
          (aset array-v (dec ny) x-idx 0.0))
        (recur (inc n))))))

(defn make-plotly-quiver-annotations [{:keys [nx ny
                                              spatial-array
                                              array-u
                                              array-v] :as _params}
                                      & {:as plot-params}]
  (let [[array-x array-y] spatial-array
        !annotation (transient [])
        {:keys [step scale]} plot-params]
    (dotimes [y-idx (int (/ ny step))]
      (dotimes [x-idx (int (/ nx step))]
        (let [y-idx   (* step y-idx)
              x-idx   (* step x-idx)
              x-arrow (aget array-x x-idx)
              y-arrow (aget array-y y-idx)
              u-val   (aget array-u y-idx x-idx)
              v-val   (aget array-v y-idx x-idx)
              x-val   (+ x-arrow (* u-val scale))
              y-val   (+ y-arrow (* v-val scale))]
          (conj! !annotation {:x          x-val
                              :y          y-val
                              :ax         x-arrow
                              :ay         y-arrow
                              :xref       "x"
                              :yref       "y"
                              :axref      "x"
                              :ayref      "y"
                              :arrowhead  2
                              :arrowsize  1
                              :arrowwidth 1.5
                              :arrowcolor "blue"
                              :opacity    0.5}))))
    (persistent! !annotation)))

(def plot-params-default
  {:step   2
   :scale  0.5 ;; arrow sizes to adjust from the flow value
   :width  600
   :height 400})

(defn plotly-contour-quiver-plot [{:keys [nx ny
                                          spatial-array
                                          array-p
                                          array-u
                                          array-v] :as simulated-results}
                                  & {:as plot-params'}]
  (let [[array-x array-y] spatial-array
        contour-base-map {:x          array-x
                          :y          array-y
                          :z          array-p
                          :type       "contour"
                          :colorscale "Viridis"
                          :contours   {:coloring "fill"}
                          :opacity    0.5
                          :showscale  true}
        {:keys [title width height] :as plot-params}
        (merge plot-params-default plot-params')]
    (kind/plotly
      {:data   [(assoc contour-base-map :name "Pressure Field")
                (-> contour-base-map
                    (assoc :name "Pressure Contours")
                    (assoc-in [:contours :coloring] "line"))]
       :layout {:title       title
                :xaxis       {:title "X"}
                :yaxis       {:title "Y"}
                :width       width
                :height      height
                :annotations (make-plotly-quiver-annotations simulated-results plot-params)}})))
;;; nt = 100
;^:kindly/hide-code
;(plotly-contour-quiver-plot (cavity-flow (assoc init-params :nt 10)))

(comment
  (cavity-flow (assoc init-params :nt 18))

  (clojure.pprint/pprint array-b)
  (clojure.pprint/pprint array-u)
  (clojure.pprint/pprint array-v)
  (clojure.pprint/pprint array-p))

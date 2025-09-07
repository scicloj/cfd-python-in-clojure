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
;; The momentum equation in vector form for a velocity field $\vec{v}$ is:
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
(def dx (double (/ (- x-end x-start) (dec nx))))
(def dy (double (/ (- y-end y-start) (dec ny))))
(def spatial-init-params
  {:nx nx :x-start x-start :x-end x-end :dx dx :x-end-idx (dec nx) :x-second-end-idx (- nx 2)
   :ny ny :y-start y-start :y-end y-end :dy dy :y-end-idx (dec ny) :y-second-end-idx (- ny 2)})

(def spatial-array (two-d/create-array-2d spatial-init-params))
;;
;; The rest setup:
(def rho 1.0)
(def nu 0.1)
(def dt 0.001)

;; Doing some pre-calculation to optimize the performance
(def pre-calculation-params
  {:two-dx                (* 2.0 dx)
   :two-dy                (* 2.0 dy)
   :two-dx-dy             (* 2.0 dx dy)
   :dx-square             (pow dx 2)
   :dy-square             (pow dy 2)
   :four-dx-square        (* 4.0 (pow dx 2))
   :four-dy-square        (* 4.0 (pow dy 2))
   :dx-dy-squares         (* (pow dx 2) (pow dy 2))
   :dx-dy                 (* dx dy)
   :over-dx               (/ 1.0 dx)
   :over-dy               (/ 1.0 dy)
   :dt-over-dx            (/ dt dx)
   :dt-over-dy            (/ dt dy)
   :dt-over-dx-square     (/ dt (pow dx 2))
   :dt-over-dy-square     (/ dt (pow dy 2))
   :two-dx-dy-squares-sum (* 2.0 (+ (pow dx 2) (pow dy 2)))
   :dx-dy-squares-over-two-dx-dy-squares-sum
   (/ (* (pow dx 2) (pow dy 2)) (* 2.0 (+ (pow dx 2) (pow dy 2))))
   :dt-over-two-rho-dx    (/ dt (* 2.0 rho dx))
   :dt-over-two-rho-dy    (/ dt (* 2.0 rho dy))})

(def array-u (two-d/create-init-zeros-u (assoc spatial-init-params :d-type Double/TYPE) spatial-array))
(def array-v (two-d/clone-2d-array array-u))
(def array-p (two-d/clone-2d-array array-u))
(def array-b (two-d/clone-2d-array array-u))

;; Now we define initial parameters:
(def init-params
  (-> (merge spatial-init-params pre-calculation-params)
      (assoc :rho rho
             :nu nu
             :dt dt
             :nit nit
             :nt nt
             :spatial-array spatial-array
             :array-u array-u
             :array-v array-v
             :array-p array-p
             :array-b array-b)))
;;
;; The pressure Poisson equation that's written above can be hard to write out without typos.
;; The function `build-up-b` below represents the contents of the square brackets, so that the entirety of
;; the PPE is slightly more manageable.
;;
(defn build-up-b [{:keys [array-b array-u array-v rho dt dx dy nx ny
                          x-end-idx y-end-idx
                          x-second-end-idx y-second-end-idx
                          two-dx
                          two-dy
                          two-dx-dy
                          four-dx-square
                          four-dy-square] :as _params}]
  (dotimes [y-idx y-second-end-idx]
    (let [y-prev-idx y-idx
          y-idx      (inc y-idx)
          y-next-idx (inc y-idx)

          u-prev-x   (aget array-u y-prev-idx)
          u-x        (aget array-u y-idx)
          u-next-x   (aget array-u y-next-idx)
          v-prev-x   (aget array-v y-prev-idx)
          v-x        (aget array-v y-idx)
          v-next-x   (aget array-v y-next-idx)]
      (dotimes [x-idx x-second-end-idx]
        (let [;; index definitions
              x-prev-idx x-idx
              x-idx      (inc x-idx)
              x-next-idx (inc x-idx)

              ;; extracting array values
              u-i-j-1    (aget u-prev-x x-idx)
              u-i-j+1    (aget u-next-x x-idx)
              u-j-i-1    (aget u-x x-prev-idx)
              u-j-i+1    (aget u-x x-next-idx)

              v-i-j-1    (aget v-prev-x x-idx)
              v-i-j+1    (aget v-next-x x-idx)
              v-j-i-1    (aget v-x x-prev-idx)
              v-j-i+1    (aget v-x x-next-idx)

              ;; pre-calculations
              u-i-diff   (- u-j-i+1 u-j-i-1)
              u-j-diff   (- u-i-j+1 u-i-j-1)
              v-i-diff   (- v-j-i+1 v-j-i-1)
              v-j-diff   (- v-i-j+1 v-i-j-1)]
          (aset array-b y-idx x-idx
            (double (* rho
                       (- (/ (+ (/ u-i-diff two-dx) (/ v-j-diff two-dy)) dt)
                          (/ (* u-i-diff u-i-diff) four-dx-square)
                          (/ (* u-j-diff v-i-diff) two-dx-dy)
                          (/ (* v-j-diff v-j-diff) four-dy-square))))))))))
;;
;; The function `pressure-poisson` is also defined to help segregate the different rounds of calculations.
;; Note the presence of the pseudo-time variable `nit`.
;; This sub-iteration in the Poisson calculation helps
;; ensure a divergence-free field.
;;

(defn pressure-poisson [{:keys [array-p array-b nx ny dx dy nit
                                x-end-idx y-end-idx
                                x-second-end-idx y-second-end-idx
                                dx-square
                                dy-square
                                two-dx-dy-squares-sum
                                dx-dy-squares-over-two-dx-dy-squares-sum]
                         :as   _params}]
  (dotimes [_ nit]
    (let [pn (two-d/clone-2d-array array-p)]
      (dotimes [y-idx y-second-end-idx]
        (let [y-prev-idx y-idx
              y-idx      (inc y-idx)
              y-next-idx (inc y-idx)

              pn-x       (aget pn y-idx)
              pn-prev-x  (aget pn y-prev-idx)
              pn-next-x  (aget pn y-next-idx)
              array-b-x  (aget array-b y-idx)]
          (dotimes [x-idx x-second-end-idx]
            (let [x-prev-idx   x-idx
                  x-idx        (inc x-idx)
                  x-next-idx   (inc x-idx)

                  p-i-j-1      (aget pn-prev-x x-idx)
                  p-i-j+1      (aget pn-next-x x-idx)
                  p-j-i-1      (aget pn-x x-prev-idx)
                  p-j-i+1      (aget pn-x x-next-idx)
                  p-i-diff-sum (+ p-j-i+1 p-j-i-1)
                  p-j-diff-sum (+ p-i-j+1 p-i-j-1)]
              (aset array-p y-idx x-idx
                (double (- (/ (+ (* p-i-diff-sum dy-square)
                                 (* p-j-diff-sum dx-square))
                              two-dx-dy-squares-sum)
                           (* dx-dy-squares-over-two-dx-dy-squares-sum
                              (aget array-b-x x-idx)))))))))

      ;; boundary conditions
      ;; dp/dx = 0 at x = 2 & dp/dx = 0 at x = 0
      (dotimes [y-idx ny]
        (aset array-p y-idx x-end-idx (aget array-p y-idx x-second-end-idx))
        (aset array-p y-idx 0 (aget array-p y-idx 1)))

      ;; dp/dy = 0 at y = 0 & p = 0 at y = 2
      (dotimes [x-idx nx]
        (aset array-p 0 x-idx (aget array-p 1 x-idx))
        (aset array-p y-end-idx x-idx 0.0)))))

;;
;; Finally, the rest of the cavity flow equations are wrapped inside the function `cavity-flow`,
;; allowing us to easily plot the results of the cavity flow solver for different lengths of time.
;;
(def !nt (atom 0))

(defn cavity-flow [{:keys [nt array-u array-v array-p dt nx ny dx dy rho nu
                           x-end-idx y-end-idx
                           x-second-end-idx y-second-end-idx
                           dt-over-dx
                           dt-over-dy
                           dt-over-dx-square
                           dt-over-dy-square
                           two-dx-dy-squares-sum
                           dt-over-two-rho-dx
                           dt-over-two-rho-dy] :as params}]
  (while (< @!nt nt)
   (let [un (two-d/clone-2d-array array-u)
         vn (two-d/clone-2d-array array-v)]
     (build-up-b params)
     (pressure-poisson params)

     (dotimes [y-idx y-second-end-idx]
       (let [y-prev-idx y-idx
             y-idx      (inc y-idx)
             y-next-idx (inc y-idx)

             un-prev-x  (aget un y-prev-idx)
             un-x       (aget un y-idx)
             un-next-x  (aget un y-next-idx)

             vn-prev-x  (aget vn y-prev-idx)
             vn-x       (aget vn y-idx)
             vn-next-x  (aget vn y-next-idx)

             p-prev-x   (aget array-p y-prev-idx)
             p-x        (aget array-p y-idx)
             p-next-x   (aget array-p y-next-idx)]
        (dotimes [x-idx x-second-end-idx]
          (let [x-prev-idx    x-idx
                x-idx         (inc x-idx)
                x-next-idx    (inc x-idx)

                u-i-j         (aget un-x x-idx)
                u-i-j-1       (aget un-prev-x x-idx)
                u-i-j+1       (aget un-next-x x-idx)
                u-j-i-1       (aget un-x x-prev-idx)
                u-j-i+1       (aget un-x x-next-idx)
                neg-two-u-i-j (* -2.0 u-i-j)

                v-i-j         (aget vn-x x-idx)
                v-i-j-1       (aget vn-prev-x x-idx)
                v-i-j+1       (aget vn-next-x x-idx)
                v-j-i-1       (aget vn-x x-prev-idx)
                v-j-i+1       (aget vn-x x-next-idx)
                neg-two-v-i-j (* -2.0 v-i-j)

                p-j-i+1       (aget p-x x-next-idx)
                p-j-i-1       (aget p-x x-prev-idx)
                p-i-j+1       (aget p-next-x x-idx)
                p-i-j-1       (aget p-prev-x x-idx)]
            (aset array-u y-idx x-idx
              (double (- u-i-j
                         (* u-i-j dt-over-dx (- u-i-j u-j-i-1))
                         (* v-i-j dt-over-dy (- u-i-j u-i-j-1))
                         (* dt-over-two-rho-dx (- p-j-i+1 p-j-i-1))
                         (* (- nu)
                            (+ (* dt-over-dx-square
                                  (+ u-j-i+1 neg-two-u-i-j u-j-i-1))
                               (* dt-over-dy-square
                                  (+ u-i-j+1 neg-two-u-i-j u-i-j-1)))))))

            (aset array-v y-idx x-idx
              (double (- v-i-j
                         (* u-i-j dt-over-dx (- v-i-j v-j-i-1))
                         (* v-i-j dt-over-dy (- v-i-j v-i-j-1))
                         (* dt-over-two-rho-dy (- p-i-j+1 p-i-j-1))
                         (* (- nu)
                            (+ (* dt-over-dx-square
                                  (+ v-j-i+1 neg-two-v-i-j v-j-i-1))
                               (* dt-over-dy-square
                                  (+ v-i-j+1 neg-two-v-i-j v-i-j-1)))))))))))

     ;; boundary conditions
     (dotimes [y-idx ny]
       (aset array-u y-idx 0 0.0)
       (aset array-u y-idx x-end-idx 0.0)
       (aset array-v y-idx 0 0.0)
       (aset array-v y-idx x-end-idx 0.0))
     (dotimes [x-idx nx]
       (aset array-u 0 x-idx 0.0)
       (aset array-u y-end-idx x-idx 1.0) ;; set velocity on cavity lid equal to 1
       (aset array-v 0 x-idx 0.0)
       (aset array-v y-end-idx x-idx 0.0))
     (swap! !nt inc))))
;;
;; Run simulation with nt = 100
(cavity-flow (assoc init-params :nt 100))
;;
^:kindly/hide-code
(two-d/plotly-contour-quiver-plot init-params)

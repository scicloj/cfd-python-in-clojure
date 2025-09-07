^:kindly/hide-code
(ns steps.step-12
  (:require
   [cfd.two-d :as two-d]
   [fastmath.core :as fm :refer [pow]]
   [scicloj.kindly.v4.api :as kindly]
   [scicloj.kindly.v4.kind :as kind]
   [utils.notebook :refer [md tex]]))

;; Did you make it this far? This is the last step!
;; How long did it take you to write your own Navier-stoke solver in Clojure
;; following the interactive module? Let us know!
;;
;; ## Step 12: Channel Flow with Navier-Stokes
;;
;; The only difference between this final step and Step 11 is that we are
;; going to add a source term to the $u$-momentum equation, to mimic the effect
;; of a pressure-driven channel flow. Here are our modified Navier-Stokes equations:
;;
(tex "\\frac{\\partial u}{\\partial t}+u\\frac{\\partial u}{\\partial x}+v\\frac{\\partial u}{\\partial y}=-\\frac{1}{\\rho}\\frac{\\partial p}{\\partial x}+\\nu\\left(\\frac{\\partial^2 u}{\\partial x^2}+\\frac{\\partial^2 u}{\\partial y^2}\\right)+F")
(tex "\\frac{\\partial v}{\\partial t}+u\\frac{\\partial v}{\\partial x}+v\\frac{\\partial v}{\\partial y}=-\\frac{1}{\\rho}\\frac{\\partial p}{\\partial y}+\\nu\\left(\\frac{\\partial^2 v}{\\partial x^2}+\\frac{\\partial^2 v}{\\partial y^2}\\right)")
(tex "\\frac{\\partial^2 p}{\\partial x^2}+\\frac{\\partial^2 p}{\\partial y^2}=-\\rho\\left(\\frac{\\partial u}{\\partial x}\\frac{\\partial u}{\\partial x}+2\\frac{\\partial u}{\\partial y}\\frac{\\partial v}{\\partial x}+\\frac{\\partial v}{\\partial y}\\frac{\\partial v}{\\partial y}\\right)")
;;
;; ### Discretized equations
;;
;; With patience and care, we write the discretized form of the equations.
;; It is highly recommended that you write these in your own hand, mentally
;; following each term as you write it.
;;
;; The $u$-momentum equation:
(tex "\\begin{split}\n& \\frac{u_{i,j}^{n+1}-u_{i,j}^{n}}{\\Delta t}+u_{i,j}^{n}\\frac{u_{i,j}^{n}-u_{i-1,j}^{n}}{\\Delta x}+v_{i,j}^{n}\\frac{u_{i,j}^{n}-u_{i,j-1}^{n}}{\\Delta y} = \\\\\n& \\qquad -\\frac{1}{\\rho}\\frac{p_{i+1,j}^{n}-p_{i-1,j}^{n}}{2\\Delta x} \\\\\n& \\qquad +\\nu\\left(\\frac{u_{i+1,j}^{n}-2u_{i,j}^{n}+u_{i-1,j}^{n}}{\\Delta x^2}+\\frac{u_{i,j+1}^{n}-2u_{i,j}^{n}+u_{i,j-1}^{n}}{\\Delta y^2}\\right)+F_{i,j}\n\\end{split}")
;;
;; The $v$-momentum equation:
(tex "\\begin{split}\n& \\frac{v_{i,j}^{n+1}-v_{i,j}^{n}}{\\Delta t}+u_{i,j}^{n}\\frac{v_{i,j}^{n}-v_{i-1,j}^{n}}{\\Delta x}+v_{i,j}^{n}\\frac{v_{i,j}^{n}-v_{i,j-1}^{n}}{\\Delta y} = \\\\\n& \\qquad -\\frac{1}{\\rho}\\frac{p_{i,j+1}^{n}-p_{i,j-1}^{n}}{2\\Delta y} \\\\\n& \\qquad +\\nu\\left(\\frac{v_{i+1,j}^{n}-2v_{i,j}^{n}+v_{i-1,j}^{n}}{\\Delta x^2}+\\frac{v_{i,j+1}^{n}-2v_{i,j}^{n}+v_{i,j-1}^{n}}{\\Delta y^2}\\right)\n\\end{split}")
;;
;; And the pressure equation:
(tex "\\begin{split}\n& \\frac{p_{i+1,j}^{n}-2p_{i,j}^{n}+p_{i-1,j}^{n}}{\\Delta x^2} + \\frac{p_{i,j+1}^{n}-2p_{i,j}^{n}+p_{i,j-1}^{n}}{\\Delta y^2} = \\\\\n& \\qquad \\rho\\left[\\frac{1}{\\Delta t}\\left(\\frac{u_{i+1,j}-u_{i-1,j}}{2\\Delta x}+\\frac{v_{i,j+1}-v_{i,j-1}}{2\\Delta y}\\right) - \\frac{u_{i+1,j}-u_{i-1,j}}{2\\Delta x}\\frac{u_{i+1,j}-u_{i-1,j}}{2\\Delta x} - 2\\frac{u_{i,j+1}-u_{i,j-1}}{2\\Delta y}\\frac{v_{i+1,j}-v_{i-1,j}}{2\\Delta x} - \\frac{v_{i,j+1}-v_{i,j-1}}{2\\Delta y}\\frac{v_{i,j+1}-v_{i,j-1}}{2\\Delta y}\\right]\n\\end{split}")
;;
;; As always, we need to re-arrange these equations to the form we need in the code
;; to make the iterations proceed.
;;
;; For the $u$ and $v$ momentum equations, we isolate the velocity at the time step `n+1`:
(tex "\\begin{split}\nu_{i,j}^{n+1} = u_{i,j}^{n} & - u_{i,j}^{n} \\frac{\\Delta t}{\\Delta x} \\left(u_{i,j}^{n}-u_{i-1,j}^{n}\\right) - v_{i,j}^{n} \\frac{\\Delta t}{\\Delta y} \\left(u_{i,j}^{n}-u_{i,j-1}^{n}\\right) \\\\\n& - \\frac{\\Delta t}{\\rho 2\\Delta x} \\left(p_{i+1,j}^{n}-p_{i-1,j}^{n}\\right) \\\\\n& + \\nu\\left[\\frac{\\Delta t}{\\Delta x^2} \\left(u_{i+1,j}^{n}-2u_{i,j}^{n}+u_{i-1,j}^{n}\\right) + \\frac{\\Delta t}{\\Delta y^2} \\left(u_{i,j+1}^{n}-2u_{i,j}^{n}+u_{i,j-1}^{n}\\right)\\right] \\\\\n& + \\Delta t F\n\\end{split}")
(tex "\\begin{split}\nv_{i,j}^{n+1} = v_{i,j}^{n} & - u_{i,j}^{n} \\frac{\\Delta t}{\\Delta x} \\left(v_{i,j}^{n}-v_{i-1,j}^{n}\\right) - v_{i,j}^{n} \\frac{\\Delta t}{\\Delta y} \\left(v_{i,j}^{n}-v_{i,j-1}^{n}\\right) \\\\\n& - \\frac{\\Delta t}{\\rho 2\\Delta y} \\left(p_{i,j+1}^{n}-p_{i,j-1}^{n}\\right) \\\\\n& + \\nu\\left[\\frac{\\Delta t}{\\Delta x^2} \\left(v_{i+1,j}^{n}-2v_{i,j}^{n}+v_{i-1,j}^{n}\\right) + \\frac{\\Delta t}{\\Delta y^2} \\left(v_{i,j+1}^{n}-2v_{i,j}^{n}+v_{i,j-1}^{n}\\right)\\right]\n\\end{split}")
;;
;; And for the pressure equation, we isolate the term $p_{i,j}^n$ to iterate in pseudo-time:
(tex "\\begin{split}\np_{i,j}^{n} = & \\frac{\\left(p_{i+1,j}^{n}+p_{i-1,j}^{n}\\right) \\Delta y^2 + \\left(p_{i,j+1}^{n}+p_{i,j-1}^{n}\\right) \\Delta x^2}{2(\\Delta x^2+\\Delta y^2)} \\\\\n& -\\frac{\\rho\\Delta x^2\\Delta y^2}{2\\left(\\Delta x^2+\\Delta y^2\\right)} \\\\\n& \\times \\left[\\frac{1}{\\Delta t} \\left(\\frac{u_{i+1,j}-u_{i-1,j}}{2\\Delta x} + \\frac{v_{i,j+1}-v_{i,j-1}}{2\\Delta y}\\right) - \\frac{u_{i+1,j}-u_{i-1,j}}{2\\Delta x}\\frac{u_{i+1,j}-u_{i-1,j}}{2\\Delta x} - 2\\frac{u_{i,j+1}-u_{i,j-1}}{2\\Delta y}\\frac{v_{i+1,j}-v_{i-1,j}}{2\\Delta x} - \\frac{v_{i,j+1}-v_{i,j-1}}{2\\Delta y}\\frac{v_{i,j+1}-v_{i,j-1}}{2\\Delta y}\\right]\n\\end{split}")
;;
;; The initial condition is $u, v, p = 0$ everywhere, and at the boundary conditions are:
;; $u, v, p$ are periodic on $x = 0, 2$ <br>
;; $u, v = 0$ at $y = 0, 2$ <br>
;; $\frac{\partial p}{\partial y} = 0$ at $y = 0, 2$ <br>
;; $F = 1$ everywhere
;;
;; In step 11, we isolated a portion of our transposed equation to make it easier to parse
;; and we're going to do the same thing here. One thing to note is that we have periodic
;; boundary conditions throughout this grid, so we need to explicitly calculate the values
;; at the leading and trailing edge of our `u` vector.
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
    (let [y-prev-idx                 y-idx
          y-idx                      (inc y-idx)
          y-next-idx                 (inc y-idx)

          u-prev-x                   (aget array-u y-prev-idx)
          u-x                        (aget array-u y-idx)
          u-next-x                   (aget array-u y-next-idx)
          v-prev-x                   (aget array-v y-prev-idx)
          v-x                        (aget array-v y-idx)
          v-next-x                   (aget array-v y-next-idx)

          u-x-start                  (aget u-x 0)
          u-x-2nd                    (aget u-x 1)
          u-x-2nd-end                (aget u-x x-second-end-idx)
          u-x-end                    (aget u-x x-end-idx)
          u-next-x-start             (aget u-next-x 0)
          u-prev-x-start             (aget u-prev-x 0)
          u-next-x-end               (aget u-next-x x-end-idx)
          u-prev-x-end               (aget u-prev-x x-end-idx)

          v-x-start                  (aget v-x 0)
          v-x-2nd                    (aget v-x 1)
          v-x-2nd-end                (aget v-x x-second-end-idx)
          v-x-end                    (aget v-x x-end-idx)
          v-next-x-start             (aget v-next-x 0)
          v-prev-x-start             (aget v-prev-x 0)
          v-next-x-end               (aget v-next-x x-end-idx)
          v-prev-x-end               (aget v-prev-x x-end-idx)

          u-x-start-end-diff         (- u-x-start u-x-2nd-end)
          u-x-start-end-diff-sqr     (* u-x-start-end-diff u-x-start-end-diff)

          v-next-prev-end-diff       (- v-next-x-end v-prev-x-end)
          v-next-prev-end-diff-sqr   (* v-next-prev-end-diff v-next-prev-end-diff)

          u-x-2nd-end-diff           (- u-x-2nd u-x-end)
          u-x-2nd-end-diff-sqr       (* u-x-2nd-end-diff u-x-2nd-end-diff)

          v-next-prev-start-diff     (- v-next-x-start v-prev-x-start)
          v-next-prev-start-diff-sqr (* v-next-prev-start-diff v-next-prev-start-diff)]
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
                          (/ (* v-j-diff v-j-diff) four-dy-square)))))))

      ;; Periodic BC Pressure @ x = 2
      (aset array-b y-idx x-end-idx
        (double (* rho
                   (- (/ (+ (/ u-x-start-end-diff two-dx)
                            (/ v-next-prev-end-diff two-dy))
                         dt)
                      (/ u-x-start-end-diff-sqr four-dx-square)
                      (/ (* (- u-next-x-end u-prev-x-end) (- v-x-start v-x-2nd-end)) two-dx-dy)
                      (/ v-next-prev-end-diff-sqr four-dy-square)))))

      ;; Periodic BC Pressure @ x = 0
      (aset array-b y-idx 0
        (double (* rho
                   (- (/ (+ (/ u-x-2nd-end-diff two-dx)
                            (/ v-next-prev-start-diff two-dy))
                         dt)
                      (/ u-x-2nd-end-diff-sqr four-dx-square)
                      (/ (* (- u-next-x-start u-prev-x-start) (- v-x-2nd v-x-end)) two-dx-dy)
                      (/ v-next-prev-start-diff-sqr four-dy-square))))))))
;;
;; We'll also define a Pressure Poisson iterative function, again like we did in Step 11.
;; Once more, note that we have to include the periodic boundary conditions at the leading and
;; trailing edge.
;; We also have to specify the boundary conditions at the top and bottom of our grid.
;;
(defn pressure-poisson-periodic [{:keys [array-p array-b nx ny dx dy nit
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
                              (aget array-b-x x-idx)))))))

          ;; Periodic BC Pressure @ x = 2
          (aset array-p y-idx x-end-idx
            (double (- (/ (+ (* (+ (aget pn-x 0) (aget pn-x x-second-end-idx)) dy-square)
                             (* (+ (aget pn-next-x x-end-idx) (aget pn-prev-x x-end-idx)) dx-square))
                          two-dx-dy-squares-sum)
                       (* dx-dy-squares-over-two-dx-dy-squares-sum
                          (aget array-b-x x-end-idx)))))

          ;; Periodic BC Pressure @ x = 0
          (aset array-p y-idx 0
            (double (- (/ (+ (* (+ (aget pn-x 1) (aget pn-x x-end-idx)) dy-square)
                             (* (+ (aget pn-next-x 0) (aget pn-prev-x 0)) dx-square))
                          two-dx-dy-squares-sum)
                       (* dx-dy-squares-over-two-dx-dy-squares-sum
                          (aget array-b-x 0)))))))

      ;; Wall boundary conditions, pressure
      (aset array-p y-end-idx (aget array-p y-second-end-idx))
      (aset array-p 0 (aget array-p 1)))))
;;
;; Now we have our familiar list of variables and initial conditions to declare before we start.
;;
;; ### Variable Declarations
;;
;; **Spatial setup:**
(def nx 41)
(def ny 41)
(def nt 10)
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
;;
(def spatial-array (two-d/create-array-2d spatial-init-params))
;;
;; **The rest setup:**
(def rho 1.0)
(def nu 1e-1)
(def F 1.0)
(def dt 1e-2)
;;
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
   :dt-over-two-rho-dy    (/ dt (* 2.0 rho dy))
   :dt-F                  (* dt F)})

(def array-u (two-d/create-init-zeros-u (assoc spatial-init-params :d-type Double/TYPE) spatial-array))
(def array-v (two-d/clone-2d-array array-u))
(def array-p (two-d/clone-2d-array array-u))
(def array-b (two-d/clone-2d-array array-u))
;;
;; Now we define initial parameters:
(def init-params
  (-> (merge spatial-init-params pre-calculation-params)
      (assoc :rho rho
             :nu nu
             :dt dt
             :nit nit
             :nt nt
             :F F
             :spatial-array spatial-array
             :array-u array-u
             :array-v array-v
             :array-p array-p
             :array-b array-b)))
;;
;; For the meat of our computation, we're going to reach back to a trick we used in Step 9
;; for Laplace's Equation. We're interested in what our grid will look like once we've
;; reached a near-steady state. We can either specify a number of timesteps `nt` and
;; increment it until we're satisfied with the results, or we can tell our code to run until the
;; difference between two consecutive iterations is very small.
;;
;; We also have to manage **8** separate boundary conditions for each iteration. The code
;; below writes each of them out explicitly. If you're interested in a challenge, you can try
;; to write a function which can handle some or all these boundary conditions.
;; At the moment, I am just trying to make the simulation running, so - todo: refactor it to be better
;;
(def default-u-diff 1.0)
(def !step-count (atom 0))
(def !u-diff (atom default-u-diff))
;;
(defn cavity-flow-periodic [{:keys [nt array-u array-v array-p dt nx ny dx dy rho nu
                                    x-end-idx y-end-idx
                                    x-second-end-idx y-second-end-idx
                                    dt-over-dx
                                    dt-over-dy
                                    dt-over-dx-square
                                    dt-over-dy-square
                                    two-dx-dy-squares-sum
                                    dt-over-two-rho-dx
                                    dt-over-two-rho-dy
                                    dt-F]
                             :as   params}]
  (let [neg-dt-F (- dt-F)]
    (while (> @!u-diff 1e-3)
      (let [un (two-d/clone-2d-array array-u)
            vn (two-d/clone-2d-array array-v)]
        (build-up-b params)
        (pressure-poisson-periodic params)
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
                                      (+ u-i-j+1 neg-two-u-i-j u-i-j-1))))
                             neg-dt-F)))

                (aset array-v y-idx x-idx
                  (double (- v-i-j
                             (* u-i-j dt-over-dx (- v-i-j v-j-i-1))
                             (* v-i-j dt-over-dy (- v-i-j v-i-j-1))
                             (* dt-over-two-rho-dy (- p-i-j+1 p-i-j-1))
                             (* (- nu)
                                (+ (* dt-over-dx-square
                                      (+ v-j-i+1 neg-two-v-i-j v-j-i-1))
                                   (* dt-over-dy-square
                                      (+ v-i-j+1 neg-two-v-i-j v-i-j-1)))))))))

            ;; Boundary conditions
            (let [un-x-end        (aget un-x x-end-idx)
                  un-x-2nd-end    (aget un-x x-second-end-idx)
                  un-x-start      (aget un-x 0)
                  un-x-2nd        (aget un-x 1)
                  un-prev-x-start (aget un-prev-x 0)
                  un-prev-x-end   (aget un-prev-x x-end-idx)
                  un-next-x-end   (aget un-next-x x-end-idx)
                  un-next-x-start (aget un-next-x 0)

                  vn-x-end        (aget vn-x x-end-idx)
                  vn-x-2nd-end    (aget vn-x x-second-end-idx)
                  vn-x-start      (aget vn-x 0)
                  vn-x-2nd        (aget vn-x 1)
                  vn-prev-x-start (aget vn-prev-x 0)
                  vn-prev-x-end   (aget vn-prev-x x-end-idx)
                  vn-next-x-end   (aget vn-next-x x-end-idx)
                  vn-next-x-start (aget vn-next-x 0)

                  p-x-start       (aget p-x 0)
                  p-x-2nd         (aget p-x 1)
                  p-x-end         (aget p-x x-end-idx)
                  p-x-2nd-end     (aget p-x x-second-end-idx)
                  p-next-x-start  (aget p-next-x 0)
                  p-next-x-end    (aget p-next-x x-end-idx)
                  p-prev-x-start  (aget p-prev-x 0)
                  p-prev-x-end    (aget p-prev-x x-end-idx)]

              ;; Periodic BC u @ x = 2
              (aset array-u y-idx x-end-idx
                (double (- un-x-end
                           (* un-x-end dt-over-dx (- un-x-end un-x-2nd-end))
                           (* vn-x-end dt-over-dy (- un-x-end un-prev-x-end))
                           (* dt-over-two-rho-dx (- p-x-start p-x-2nd-end))
                           (* (- nu)
                              (+ (* dt-over-dx-square
                                    (+ un-x-start (* -2.0 un-x-end) un-x-2nd-end))
                                 (* dt-over-dy-square
                                    (+ un-next-x-end (* -2.0 un-x-end) un-prev-x-end))))
                           neg-dt-F)))

              ;; Periodic BC u @ x = 0
              (aset array-u y-idx 0
                (double (- un-x-start
                           (* un-x-start dt-over-dx (- un-x-start un-x-end))
                           (* vn-x-start dt-over-dy (- un-x-start un-prev-x-start))
                           (* dt-over-two-rho-dx (- p-x-2nd p-x-end))
                           (* (- nu)
                              (+ (* dt-over-dx-square
                                    (+ un-x-2nd (* -2.0 un-x-start) un-x-end))
                                 (* dt-over-dy-square
                                    (+ un-next-x-start (* -2.0 un-x-start) un-prev-x-start))))
                           neg-dt-F)))

              ;; Periodic BC v @ x = 2
              (aset array-v y-idx x-end-idx
                (double (- vn-x-end
                           (* un-x-end dt-over-dx (- vn-x-end vn-x-2nd-end))
                           (* vn-x-end dt-over-dy (- vn-x-end vn-prev-x-end))
                           (* dt-over-two-rho-dy (- p-next-x-end p-prev-x-end))
                           (* (- nu)
                              (+ (* dt-over-dx-square
                                    (+ vn-x-start (* -2.0 vn-x-end) vn-x-2nd-end))
                                 (* dt-over-dy-square
                                    (+ vn-next-x-end (* -2.0 vn-x-end) vn-prev-x-end)))))))

              ;; Periodic BC v @ x = 0
              (aset array-v y-idx 0
                (double (- vn-x-start
                           (* un-x-start dt-over-dx (- vn-x-start vn-x-end))
                           (* vn-x-start dt-over-dy (- vn-x-start vn-prev-x-start))
                           (* dt-over-two-rho-dy (- p-next-x-start p-prev-x-start))
                           (* (- nu)
                              (+ (* dt-over-dx-square
                                    (+ vn-x-2nd (* -2.0 vn-x-start) vn-x-end))
                                 (* dt-over-dy-square
                                    (+ vn-next-x-start (* -2.0 vn-x-start) vn-prev-x-start))))))))))

        ;; Wall boundary conditions: u,v = 0 @ y = 0, 2
        (dotimes [x-idx nx]
          (aset array-u 0 x-idx 0.0)
          (aset array-u y-end-idx x-idx 0.0)
          (aset array-v 0 x-idx 0.0)
          (aset array-v y-end-idx x-idx 0.0))

        (let [sum-u  (reduce + 0.0 (mapcat identity array-u))
              sum-un (reduce + 0.0 (mapcat identity un))
              u-diff (if (zero? sum-u) 0.0 (/ (- sum-u sum-un) sum-u))]
          (reset! !u-diff u-diff))
        (swap! !step-count inc)))))
;;
;; Run the simulation
(cavity-flow-periodic init-params)
;;
;; You can see that we've also included `!step-count` to see how many iterations our loop
;; went through before our stop condition was met.
;;
^:kindly/hide-code
(str @!step-count)
;;
;; If you want to see how the number of iteration increases as our `!u-diff` condition gets smaller,
;; try defining a function to perform the `while` loop written above that takes an input `!`u-diff`
;; and outputs the number of iterations that the function runs.
;;
;; For now, let's look at our results. We've used the quiver functions to look at the cavity flow results
;; and it works well for channel flow too.
(two-d/plotly-quiver-plot init-params :step 3 :scale 0.1)
;;
;; `:step 3` are useful when dealing with large amount of data that you want to visualize.
;; The one used above tells `plotly` to only plot every 3rd data point. If we leave it out, you can see
;; that the result can appear a little crowded since I set the default value `2` in the plotting function.
;;
;; ### Learn more
;;
;; **What is the meaning of the $F$ term?**
;;
;; Step 12 is an exercise demonstrating the problem of flow in a channel or pipe. If you recall from fluid mechanics
;; class, a specified pressure gradient is what drives Poisseulle flow.
;;
;; Recall the $x$-momentum equation:
(tex "\\frac{\\partial u}{\\partial t}+u \\cdot \\nabla u = -\\frac{\\partial p}{\\partial x}+\\nu \\nabla^2 u")
;;
;; What we actually do in Step 12 is split the pressure into steady and unsteady components $p=P+p'$.
;; The applied steady pressure gradient is the constant $-\frac{\partial P}{\partial x}=F$(interpreted as a source term),
;; and the unsteady component is $\frac{\partial p'}{\partial x}$. So the pressure that we solve for in Step 12 is
;; actually $p'$, which for a steady flow is in fact equal to zero everywhere.
;;
;; **Why did we do this?**
;;
;; Note that we used periodic boundary conditions for this flow. For a flow with a constant pressure gradient,
;; the value of pressure on the left edge of the domain must be different from the pressure at the right edge.
;; So we cannot apply periodic boundary conditions on the pressure directly.
;; It is easier to fix the gradient and solve the perturbations in pressure.
;;
;; **Shouldn't we always expect a uniform/constant $p'$ then?**
;;
;; That's true only in the case of steady laminar flows. At high Reynolds numbers, flows in channels can become
;; turbulent, and we will see unsteady fluctuations in the pressure, which will result in non-zero values for $p'$.
;;
;; In step 12, note that the pressure field itself is not constant, but it's the pressure perturbation field is.
;; The pressure field varies linearly along the channel with slope equal to the pressure gradient. Also, for
;; incompressible flows, the absolute value of the pressure is inconsequential.
;;
;; **And explore more CFD materials online**
;;
;; The interactive module **12 steps to Navier-Stokes** is one of several components of the Computational Fluid Dynamics class
;; taught by Prof. Lorena A. Barba in Boston University between 2009 and 2013.
;;
;; For a sample of what the other components of this class are, you can explore the **Resources** section of the Spring 2013
;; version of [the course's Piazza site](https://piazza.com/bu/spring2013/me702/resources).
;;

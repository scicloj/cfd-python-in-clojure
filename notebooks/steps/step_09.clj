^:kindly/hide-code
(ns steps.step-09
  (:require
   [cfd.two-d :as two-d]
   [fastmath.core :as fm :refer [pow]]
   [scicloj.kindly.v4.api :as kindly]
   [scicloj.kindly.v4.kind :as kind]
   [utils.notebook :refer [md tex]]))

;; In the previous step, we solved the 2D Burgers' equation: an important equation in the study of fluid
;; mechanics because it contains the full convective nonlinearity of the flow equations. With that exercise,
;; we also build the experience to incrementally code a Navier-Stokes solver.
;;
;; In the next two steps, we will solve Laplace and then Poisson equation. We will then put it all together!
;;
;; ## Step 9: 2-D Laplace Equation
;;
;; Here is [Laplace's equation](https://en.wikipedia.org/wiki/Laplace%27s_equation) in 2D:
(tex "\\frac{\\partial ^2 p}{\\partial x^2} + \\frac{\\partial ^2 p}{\\partial y^2} = 0")
;; We know how to discretize a 2nd order derivative. But this about this for a minute - Laplace's
;; equation has the features typical of diffusion phenomena. For this reason, it has to be discretized with
;; [central differences](https://en.wikipedia.org/wiki/Central_differencing_scheme), so that the
;; discretization is consistent with the physics we want to simulate.
;;
;; The discretized equation is:
(tex "\\frac{p_{i+1, j}^n - 2p_{i,j}^n + p_{i-1,j}^n}{\\Delta x^2} + \\frac{p_{i,j+1}^n - 2p_{i,j}^n + p_{i, j-1}^n}{\\Delta y^2} = 0")
;; Notice that the Laplace Equation does not have a time dependence - there is no $p^{n+1}$.
;; Instead of tracking a wave through time (like in the previous steps), the Laplace equation calculates
;; the equilibrium state of a system under the supplied boundary conditions.
;;
;; If you have taken coursework in Heat Transfer, you will recognize the Laplace Equation
;; as the steady-state heat equation.
;;
;; Instead of calculating where the system will be at some time $t$, we will iteratively solve for
;; $p_{i,j}^n$ until it meets a condition that we specify. The system will reach equilibrium only
;; as the number of iterations tends to $\infty$, but we can approximate the equilibrium state by
;; iterating until the change between one iteration and the next is very small.
;;
;; Let's rearrange the discretized equation, solving for $p_{i,j}^n$:
(tex "p_{i,j}^n = \\frac{\\Delta y^2(p_{i+1,j}^n+p_{i-1,j}^n)+\\Delta x^2(p_{i,j+1}^n + p_{i,j-1}^n)}{2(\\Delta x^2 + \\Delta y^2)}")
;; Using second order central-difference schemes in both directions is the most widely applied method
;; for the Laplace operator. It is also known as the **five-point difference operator**, alluding to its stencil.
;;
;; We are going to solve Laplace's equation numerically by assuming an initial state of $p = 0$ everywhere.
;; Then we add boundary conditions as follows:
;;
;; $p = 0$ at $x = 0$
;; $p = y$ at $x = 2$
;; $\frac{\\Delta p}{\\Delta y} = 0$ at $y = 0, 1$
;;
;; Under these conditions, there is an analytical solution for Laplace's equation:
(tex "p(x,y)=\\frac{x}{4}-4\\sum_{n=1,odd}^{\\infty}\\frac{1}{(n\\pi)^2\\sinh2n\\pi}\\sinh n\\pi x\\cos n\\pi y")
;;
;; Define the initial parameters
(def nx 31)
(def ny 31)
(def spatial-init-param
  {:nx nx :x-start 0 :x-end 2 :dx (double (/ 2 (- nx -1)))
   :ny ny :y-start 0 :y-end 1 :dy (double (/ 2 (- ny -1)))})
;;
;; Have the discretized 2D spatial array ready:
(def spatial-array (two-d/create-array-2d spatial-init-param))

(def array-p (two-d/create-init-u
               (assoc spatial-init-param
                 :d-type Double/TYPE
                 :condition-fn (fn [x-val y-val]
                                 (if (= x-val ((comp last first) spatial-array))
                                   (double y-val)
                                   (double 0))))
               spatial-array))
;;
;; ### Define the Laplace function
;;
;; We will write the Laplace function to solve for $p$ until the change in the [L1 Norm](http://en.wikipedia.org/wiki/Norm_(mathematics)#Taxicab_norm_or_Manhattan_norm)
;; of $p$ is less that a specified value.
;;

^:kindly/hide-code
(defn l1-norm [p pn]
  (let [p-flat      (mapcat identity p)
        pn-flat     (mapcat identity pn)
        numerator   (reduce + (map #(fm/abs (- (fm/abs ^Double %1) (Math/abs ^Double %2))) p-flat pn-flat))
        denominator (reduce + (map #(fm/abs %) pn-flat))]
    (/ numerator denominator)))

(def !test (atom []))
(= (mapcat identity (mapcat identity (first @!test))) (mapcat identity (second @!test)))

(defn laplace-2d [{:keys [spatial-array array-p nx ny dx dy l1norm-target]
                   :as   params}]
  (let [[_array-x array-y] spatial-array]
    (loop [l1-norm-val 1]
      (if (< l1-norm-val l1norm-target)
        array-p
        (let [pn (aclone array-p)]
          (dotimes [y-idx ny]
            (when (and (pos? y-idx) (< y-idx (dec ny)))
              (dotimes [x-idx nx]
                (when (and (pos? x-idx) (< x-idx (dec nx)))
                  (let [p-j-i+1 (aget pn y-idx (inc x-idx))
                        p-j-i-1 (aget pn y-idx (dec x-idx))
                        p-j+1-i (aget pn (inc y-idx) x-idx)
                        p-j-1-i (aget pn (dec y-idx) x-idx)
                        dx-2    (pow dx 2)
                        dy-2    (pow dy 2)]
                    (aset array-p y-idx x-idx
                      (double (/ (+ (* dy-2 (+ p-j-i+1 p-j-i-1))
                                    (* dx-2 (+ p-j+1-i p-j-1-i)))
                                 (* 2 (+ dx-2 dy-2))))))))))

          ;; boundary conditions
          ; p = 0 @ x = 0
          (dotimes [y-idx ny]
            (aset array-p y-idx 0 (double 0)))
          ;; p = y @ x = 2
          (dotimes [y-idx ny]
            (aset array-p y-idx (dec nx) (double (aget array-y y-idx))))
          ;; dp/dy = 0 @ y = 0
          (dotimes [x-idx nx]
            (aset array-p 0 x-idx (aget array-p 1 x-idx)))
          ;; dp/dy = 0 @ y = 1
          (dotimes [x-idx nx]
            (aset array-p (dec ny) x-idx (aget array-p (- ny 2) x-idx)))
          ;; calculate l1nom
          (reset! !test [array-p pn (l1-norm array-p pn)])
          (recur (l1-norm array-p pn)))))))

;; Now let's try looking at our initial condition plot.
;;
^:kindly/hide-code
(let [plottable-data (two-d/arr->plotly-plottable-data spatial-array array-p)]
  (kind/plotly
    {:data   [(merge plottable-data (assoc two-d/plotly-opts :type :scatter3d))]
     :layout {:scene {:zaxis {:range [0.0 1.0]}}}}))
;;
;; We will have init params:
(def init-params
  (merge spatial-init-param
         {:spatial-array spatial-array
          :array-p       array-p
          :l1norm-target 1e-4}))

;;
;; Then run the simulation:

(let [simulated-data (laplace-2d init-params)
      [array-x array-y] spatial-array]
  (kind/plotly
    {:colorscale :Viridis
     :data   [(merge {:x array-x
                      :y array-y
                      :z simulated-data} (assoc two-d/plotly-opts :type :surface))]
     :layout {:scene {:xaxis {:range [0.0 2.0]}
                      :yaxis {:range [0.0 1.0]}
                      :zaxis {:range [0.0 1.0]}}}}))

;; testing with surface plot

(let [simulated-data (laplace-2d init-params)
      plottable-data (two-d/arr->plotly-plottable-data spatial-array simulated-data)]
  (kind/plotly
    {:colorscale :Viridis
     :data   [(merge plottable-data two-d/plotly-opts)]
     :layout {:scene {:xaxis {:range [0.0 2.0]}
                      :yaxis {:range [0.0 1.0]}
                      :zaxis {:range [0.0 1.0]}}}}))

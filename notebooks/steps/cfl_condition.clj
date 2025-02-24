^:kindly/hide-code
(ns steps.cfl-condition
  (:require
   [cfd.one-d :as one-d]
   [scicloj.kindly.v4.api :as kindly]
   [scicloj.kindly.v4.kind :as kind]
   [utils.notebook :refer [md tex]]))

;; # Convergence and the CFL Condition

;; initial and boundary conditions:
;;
(def init-params {:x-start 0
                  :x-end   2
                  :nx      41
                  :nt      20
                  :c       1.0
                  :dt      0.025})

;; 41 points of grid and 0.025 sec of timestep.
;;
;; Experimenting increasing the size of the grid below to see what happens:
;;

;; $nx = 41$
^:kindly/hide-code
(let [params  init-params
      array-x (one-d/create-array-x params)
      array-u (one-d/create-array-u {:array-x array-x})
      u       (one-d/simulate array-u params)]
  (kind/vega-lite
    {:mark     "line"
     :width    500 :height 300
     :encoding {:x {:field "x" :type "quantitative"}
                :y {:field "y" :type "quantitative"}}
     :data     {:values (into [] (map #(hash-map :x % :y %2) array-x u))}}))

;; $nx = 61$
^:kindly/hide-code
(let [params  (-> init-params
                  (assoc :nx 61))
      array-x (one-d/create-array-x params)
      array-u (one-d/create-array-u {:array-x array-x})
      u       (one-d/simulate array-u params)]
  (kind/vega-lite
    {:mark     "line"
     :width    500 :height 300
     :encoding {:x {:field "x" :type "quantitative"}
                :y {:field "y" :type "quantitative"}}
     :data     {:values (into [] (map #(hash-map :x % :y %2) array-x u))}}))

;; $nx = 81$
^:kindly/hide-code
(let [params  (-> init-params
                  (assoc :nx 81))
      array-x (one-d/create-array-x params)
      array-u (one-d/create-array-u {:array-x array-x})
      u       (one-d/simulate array-u params)]
  (kind/vega-lite
    {:mark     "line"
     :width    500 :height 300
     :encoding {:x {:field "x" :type "quantitative"}
                :y {:field "y" :type "quantitative"}}
     :data     {:values (into [] (map #(hash-map :x % :y %2) array-x u))}}))

;; ## Reasons for the breakage
;;
;; increasing grid size means travelling distance within a $\Delta t$ becomes grater than
;; $\Delta x$, which correlates to $nx$. In order to enforce the stability, we introduce
;; **Courant number** $\sigma_{max}$. This ensures stability with given discretization params.
;;
(tex "\\sigma = \\frac{u \\Delta t}{\\Delta x} \\le \\sigma_{max}")
;;

(def init-params' {:x-start 0
                   :x-end   2
                   :nx      41
                   :nt      20
                   :c       1.0
                   :sigma   0.5})

;; $nx = 41$
^:kindly/hide-code
(let [params  init-params'
      array-x (one-d/create-array-x params)
      array-u (one-d/create-array-u {:array-x array-x})
      u       (one-d/simulate array-u params)]
  (kind/vega-lite
    {:mark     "line"
     :width    500 :height 300
     :encoding {:x {:field "x" :type "quantitative"}
                :y {:field "y" :type "quantitative"}}
     :data     {:values (into [] (map #(hash-map :x % :y %2) array-x u))}}))

;; $nx = 61$
;;
^:kindly/hide-code
(let [params  (assoc init-params' :nx 61)
      array-x (one-d/create-array-x params)
      array-u (one-d/create-array-u {:array-x array-x})
      u       (one-d/simulate array-u params)]
  (kind/vega-lite
    {:mark     "line"
     :width    500 :height 300
     :encoding {:x {:field "x" :type "quantitative"}
                :y {:field "y" :type "quantitative"}}
     :data     {:values (into [] (map #(hash-map :x % :y %2) array-x u))}}))

;; $nx = 81$
^:kindly/hide-code
(let [params  (assoc init-params' :nx 81)
      array-x (one-d/create-array-x params)
      array-u (one-d/create-array-u {:array-x array-x})
      u       (one-d/simulate array-u params)]
  (kind/vega-lite
    {:mark     "line"
     :width    500 :height 300
     :encoding {:x {:field "x" :type "quantitative"}
                :y {:field "y" :type "quantitative"}}
     :data     {:values (into [] (map #(hash-map :x % :y %2) array-x u))}}))

;; $nx = 101$
^:kindly/hide-code
(let [params  (assoc init-params' :nx 101)
      array-x (one-d/create-array-x params)
      array-u (one-d/create-array-u {:array-x array-x})
      u       (one-d/simulate array-u params)]
  (kind/vega-lite
    {:mark     "line"
     :width    500 :height 300
     :encoding {:x {:field "x" :type "quantitative"}
                :y {:field "y" :type "quantitative"}}
     :data     {:values (into [] (map #(hash-map :x % :y %2) array-x u))}}))

;; $nx = 121$
^:kindly/hide-code
(let [params  (assoc init-params' :nx 121)
      array-x (one-d/create-array-x params)
      array-u (one-d/create-array-u {:array-x array-x})
      u       (one-d/simulate array-u params)]
  (kind/vega-lite
    {:mark     "line"
     :width    500 :height 300
     :encoding {:x {:field "x" :type "quantitative"}
                :y {:field "y" :type "quantitative"}}
     :data     {:values (into [] (map #(hash-map :x % :y %2) array-x u))}}))

;; The results show with a grid size $nx$ increases, convection travels shorter distance.
;; With a given $nt = 20$ in the init param, the time windows becomes shorter as a result
;; of increasing $nx$.
;;
^:kindly/hide-code
(ns conferences.scinoj-light-1.post-talk
  (:require
   [cfd.one-d :as one-d]
   [fastmath.core :as fm :refer [PI]]
   [scicloj.kindly.v4.api :as kindly]
   [utils.notebook :refer [md tex ts-vega-lite-plot-map vega-lite-plot]]))

^:kindly/hide-code (def nx 101)
^:kindly/hide-code (def nt 500)
^:kindly/hide-code (def nu 0.05)
^:kindly/hide-code (def dx (* 2.0 PI (/ 1 (- nx 1))))
^:kindly/hide-code (def dt (* dx nu))
^:kindly/hide-code (def x-start 0)
^:kindly/hide-code (def x-end (* 1.0 PI))
^:kindly/hide-code (def init-params' {:nx nx
                                      :nt nt
                                      :dx dx
                                      :dt dt
                                      :nu nu})

^:kindly/hide-code (def init-create-x-param {:x-start x-start :x-end x-end :nx nx})
^:kindly/hide-code (def array-x (one-d/create-array-x init-create-x-param))
^:kindly/hide-code (def init-params (assoc init-params' :array-x array-x))
;;
;;
;;
;; ## 1D Shock Interaction & simulation
;;
;; ### Initial x slicing w the below initial condition:
;;
^:kindly/hide-code init-create-x-param
;;
;; ### Initial flow velocity `u` plot
;;
^:kindly/hide-code
(def init-bump-u (one-d/create-array-u {:array-x      array-x
                                        :condition-fn #(if (and (> % 1.0) (< % 1.8)) 2 0.5)}))
;;
;; Defining initial `u`: `(fn [x] (if (< x 2.0) 1 0)}))`

^:kindly/hide-code
(vega-lite-plot {:array-x array-x :array-y init-bump-u :height 260 :width 400})

;;
;; ### Linear convection simulation
;;
^:kindly/hide-code
(def linear-convection-param (assoc init-params :c 1.0))

^:kindly/hide-code (def linear-convection-arr-u (float-array init-bump-u))

^:kindly/hide-code
(def linear-convection-u (one-d/simulate-cumulate linear-convection-arr-u linear-convection-param))

^:kindly/hide-code
(vega-lite-plot {:plot-map (ts-vega-lite-plot-map [-0.1 2.1] (assoc linear-convection-param
                                                               :cum-array-y linear-convection-u))})

;;
;; ### Non-linear convection simulation
;;
^:kindly/hide-code
(def nonlinear-convection-param (assoc init-params :co-eff :nonlinear))

^:kindly/hide-code (def convection-arr-u (float-array init-bump-u))

^:kindly/hide-code
(def non-linear-convection-u (one-d/simulate-cumulate convection-arr-u nonlinear-convection-param))

^:kindly/hide-code
(vega-lite-plot {:plot-map (ts-vega-lite-plot-map [-0.1 2.1] (assoc nonlinear-convection-param
                                                               :cum-array-y non-linear-convection-u))})

;;
;; ### Diffusion simulation
;;
^:kindly/hide-code
(def diffusion-param (assoc nonlinear-convection-param :mode :diffusion))

^:kindly/hide-code (def diffusion-arr-u (float-array init-bump-u))

^:kindly/hide-code
(def diffusion-u (one-d/simulate-cumulate diffusion-arr-u diffusion-param))

^:kindly/hide-code
(vega-lite-plot {:plot-map (ts-vega-lite-plot-map [-0.1 2.1] (assoc diffusion-param
                                                               :cum-array-y diffusion-u))})


;;
;; ### 1D Burgers' equation simulation
;;
^:kindly/hide-code
(def burgers-param (assoc nonlinear-convection-param :mode :burger))

^:kindly/hide-code (def burger-arr-u (float-array init-bump-u))

^:kindly/hide-code
(def burgers-u (one-d/simulate-cumulate burger-arr-u burgers-param))

^:kindly/hide-code
(vega-lite-plot {:plot-map (ts-vega-lite-plot-map [-0.1 2.1] (assoc burgers-param
                                                               :cum-array-y burgers-u))})

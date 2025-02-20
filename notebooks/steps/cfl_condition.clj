^:kindly/hide-code
(ns steps.cfl-condition
  (:require
   [cfd.one-d :as one-d]
   [scicloj.kindly.v4.api :as kindly]
   [scicloj.kindly.v4.kind :as kind]
   [utils.notebook :refer [md tex]]))

(def init-params {:x-start 0
                  :x-end   2
                  :nx      41
                  :nt      20
                  :c       1.0
                  :dt      0.025})

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

;; $nx = 85$
^:kindly/hide-code
(let [params  (-> init-params
                  (assoc :nx 85))
      array-x (one-d/create-array-x params)
      array-u (one-d/create-array-u {:array-x array-x})
      u       (one-d/simulate array-u params)]
  (kind/vega-lite
    {:mark     "line"
     :width    500 :height 300
     :encoding {:x {:field "x" :type "quantitative"}
                :y {:field "y" :type "quantitative"}}
     :data     {:values (into [] (map #(hash-map :x % :y %2) array-x u))}}))

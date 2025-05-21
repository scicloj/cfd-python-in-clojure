(ns utils.notebook
  (:require
   [scicloj.kindly.v4.api :as kindly]
   [scicloj.kindly.v4.kind :as kind]))

(def md (comp kindly/hide-code kind/md))

(def tex (comp kindly/hide-code kind/tex))

(defn vega-lite-plot [{:keys [array-x
                              array-y
                              x-label
                              y-label
                              x-type
                              y-type
                              data
                              mark
                              width
                              height
                              plot-map
                              encoding-opts]
                       :or   {mark   "line"
                              x-type "quantitative"
                              y-type "quantitative"
                              width  500
                              height 300}
                       :as param}]
  (kind/vega-lite
    (merge
      (or plot-map
          {:mark     mark
           :width    width :height height
           :data     (or data {:values (into [] (map #(hash-map :x % :y %2) array-x array-y))})
           :encoding (merge
                       {:x {:field "x" :type x-type :title x-label}
                        :y {:field "y" :type y-type :title y-label}}
                       encoding-opts)}))))

(defn cum-vec->vega-values [{:keys [array-x cum-array-y]}]
  (apply concat
    (map-indexed
      (fn [idx array-u]
        (map #(hash-map :idx idx :x % :y %2) array-x array-u))
      cum-array-y)))

(defn cumulated-line-plot [{:keys [array-x cum-array-y partition-size]
                            :or   {partition-size 20}
                            :as   param}]
  (let [[init-arr & rest-arr] cum-array-y
        rest-arr (->> rest-arr (partition-all partition-size) (map last))
        cum-arr  (into [init-arr] rest-arr)]
    (vega-lite-plot {:data          {:values (into [] (cum-vec->vega-values (assoc param
                                                                              :cum-array-y cum-arr)))}
                     :x-label       "X"
                     :y-label       "Flow Velocity"
                     :encoding-opts {:color {:field "idx" :type "nominal" :legend nil}}
                     :height        260
                     :width         400})))

(defn ts-vega-lite-plot-map [y-domain {:keys [array-x cum-array-y nt] :as params}]
  (let [select-dt-name "selectedDt"
        initial-dt     0
        dt-step        1
        dt-field       "idx"]
    {:data      {:values (into [] (cum-vec->vega-values params))}
     :transform [{:filter (str "datum." dt-field " == " select-dt-name)}]
     :params    [{:name  select-dt-name
                  :value initial-dt
                  :bind  {:input "range"
                          :name  "Time interval idx: "
                          :min   initial-dt
                          :max   nt
                          :step  dt-step}}]
     :mark      "line"
     :encoding  {:x {:field "x"
                     :type  "quantitative"
                     :title "X"}
                 :y {:field "y"
                     :type  "quantitative"
                     :title "Flow Velocity"
                     :scale {:domain y-domain}}}
     :height    260 :width 400}))
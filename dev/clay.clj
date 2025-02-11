(ns clay
  (:require
   [scicloj.clay.v2.api :as clay]))

(def make-config
  {:base-source-path    nil #_"notebooks" ;;files watch doesn't work
   :live-reload         true
   :format              [:quarto :html]
   :base-target-path    "docs"
   :clean-up-target-dir true
   :source-path         ["notebooks/index.clj"]
   :book                {:title "CFD Python in Clojure"}})

(comment
  (clay/config)

  (clay/stop!)

  (clay/browse!)

  (clay/make! make-config))
(ns clay
  (:require
   [clojure.edn :as edn]
   [scicloj.clay.v2.api :as clay]))

(def steps-file-path "notebooks/steps.edn")

(def full-source-path
  (into ["notebooks/index.clj"]
        (->> steps-file-path
             ((comp edn/read-string slurp))
             (map :src))))

(def make-config
  (merge
    ((comp edn/read-string slurp) "clay.edn")
    {:format              [:quarto :html]
     :base-target-path    "docs"
     :clean-up-target-dir true
     :source-path         full-source-path
     :book                {:title "CFD Python in Clojure"}}))

(def talk-make-config
  {:base-source-path nil
   :live-reload      true
   :format           [:quarto :revealjs]
   :source-path      "notebooks/conferences/scinoj_light_1/siyoung_talk.clj"
   :quarto           {:format {:revealjs {:theme :serif}}}})

(comment
  (clay/config)

  (clay/stop!)

  (clay/browse!)

  ;; notebooks build
  (clay/make! make-config)

  ;; conference talk build
  (clay/make! talk-make-config))
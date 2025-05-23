(ns utils.num-clj
  "Util functions that are inspired by numpy while going through CFDPython.

  NOTE: Not a permanent ns, and functions are simplified for purposes of this repo.")

(defn linspace
  "Returns evenly spaced numbers over a specified interval.

  todo: add a param to support different data types"
  [{:keys [start stop num endpoint?]
    :or   {num       50
           endpoint? true}}]
  (let [arr  (float-array num)
        step (/ (- stop start) (cond-> num
                                 endpoint? (dec)))]
    (dotimes [i num]
      (aset arr i (float (* i step))))
    arr))

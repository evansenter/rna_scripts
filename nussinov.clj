(ns nussinov
  (:use [clojure.set :only (union difference)]))

(def min_loop_size 3)

(defn generate_table [n]
  ; http://clj-me.cgrand.net/2009/10/15/multidim-arrays/
  (map 
    (fn [i]
      (map 
        (fn [j] (if (and (>= i 1) (<= i j) (<= j (+ i min_loop_size))) 1 nil))
        (take (inc n) (iterate inc 0))))
    (take (inc n) (iterate inc 0))))
    
; (defn flush_table [table n])
  
(defn match_pairs [rna_structure i j]
  (let [padded_structure (if (= \space (first rna_structure)) rna_structure (cons \space rna_structure))]    
    (reduce 
      (fn [hash index]
        (cond
          (= \( (nth padded_structure index))
            (assoc hash index nil)
          (= \) (nth padded_structure index))
            (assoc hash (-> (select-keys hash (map first (remove last hash))) keys sort reverse first) index)
          :default
            hash))
      {} 
      (range i (inc j)))))
      
(defn closed_pairs [pairs_hash]
  (select-keys pairs_hash (map first (filter (partial every? identity) pairs_hash))))

(defn pairs_set [rna_structure i j]
  (if (or (< j i) (> i j)) 
    #{} 
    (-> (match_pairs rna_structure i j) closed_pairs set)))
      
(defn pair_distance [rna_structure i k j]
  (let [
    reference_structure (pairs_set rna_structure i j)
    comparitive_pairings (union (pairs_set rna_structure i (dec k)) (pairs_set rna_structure (inc k) (dec j)) #{[k j]})
    ] (count (union (difference reference_structure comparitive_pairings) (difference comparitive_pairings reference_structure)))))
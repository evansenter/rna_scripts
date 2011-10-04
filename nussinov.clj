; (time (pair_distance "()))())()(().).)..)).(.()))))).)(.)(..())))))())()(.)).(..)..).()..)(.)..)))(.(()....(.(...)).()...).).(.()(.((.)..()..)))()).))..))(()((((.(()))(..((..((.()((())((.(.))(()(()(..()..((()().)(.)().(..(())))(..().((()..)().))(((((((().).(..).()).(((.()())())..(.(.(...).))..)(..))).(().()..(()()((..).).)..()..((..(.)).((.(().()((()(.).))).(().(.)..(()((...)((()(().())(())).).))..)...)((.(()(.()..)(..)(())...)..)((.((((()(((.()))((()(()(()((.(.()....()..(.))..).)((()(.)()()((...(.)(()(.).(...().).()))).))(.)()(..((.))))((..))((.()).)(.)((().()..).)((.())..).(.))((())()(...))().)())(.()(.((.((.).(..()(.(..(..((((.)()()...))).)().))(..(().()).)))(.....))())().()..)()).)(.()))((().).(()(.))())((((.())((((()))()()))...()(((((.()))(.)(((.)).(....)((()((.()(((())(.(()).(..))(()((()..(.(.()))..).().)..).)(()((.).()((.())((()()).).).(.))(.(..)(..))))).()(..(.)..(..(()...)).).().().(.).(((.((()(.)).(.)..(((..)).()....)((...(((.)(.)()..(()..(...))).(....(.)))((.())()))..)().()(.()((..(()().)(..)))()" 0 500 999))

(ns nussinov
  (:use [clojure.set :only (union difference)]))

(def min_loop_size 3)

(defn generate_table [n]
  (make-array Double/TYPE (inc n) (inc n)))
  
(defn flush_table [table]
  (dotimes [i (count table)]
    (dotimes [j (count table)]
      (deep-aset doubles table i j (if (and (<= i j) (<= j (+ i min_loop_size)) (>= i 1)) 1 0)))))
  
;;  Use a list (stack) and a set after performance testing vs. Ruby version.
  
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
    
; From http://clj-me.cgrand.net/2009/10/15/multidim-arrays/
    
(defn array? [x] (-> x class .isArray))
(defn see [x] (if (array? x) (map see x) x))
    
(defmacro deep-aget
  ([hint array idx]
    `(aget ~(vary-meta array assoc :tag hint) ~idx))
  ([hint array idx & idxs]
    `(let [a# (aget ~(vary-meta array assoc :tag 'objects) ~idx)]
       (deep-aget ~hint a# ~@idxs))))
       
(defmacro deep-aset [hint array & idxsv]
  (let [hints '{doubles double ints int}
        [v idx & sxdi] (reverse idxsv)
        idxs (reverse sxdi)
        v (if-let [h (hints hint)] (list h v) v)
        nested-array (if (seq idxs)
                       `(deep-aget ~'objects ~array ~@idxs)
                        array)
        a-sym (with-meta (gensym "a") {:tag hint})]
      `(let [~a-sym ~nested-array]
         (aset ~a-sym ~idx ~v))))
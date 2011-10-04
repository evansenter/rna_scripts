; (load-file "./nussinov.clj")
; (in-ns 'nussinov)

(ns nussinov
  (:use [clojure.set :only (union difference)]))

;; --------------------------------------------------------------------
;; Constants
;; --------------------------------------------------------------------

(def temperature (+ 37 273.15))
(def boltzmann_constant 0.0019872370936902486) ; kcal / mol / K
(def base_pair_energy (Math/pow Math/E (/ 1 (* boltzmann_constant temperature))))
(def min_loop_size 3)
(def pairings { \a [\u] \u [\a \g] \g [\c \u] \c [\g] })

;; --------------------------------------------------------------------
;; Helper functions
;; --------------------------------------------------------------------

(defn array? [x] (-> x class .isArray))
(defn see [x] (if (array? x) (map see x) x))
(defn inc_range [start stop] (range start (inc stop)))
(defn symmetric_difference [set_1 set_2] (union (difference set_1 set_2) (difference set_2 set_1)))
(defn pad_string [string] (if (= \space (first string)) string (str \space string)))

;; --------------------------------------------------------------------
; From http://clj-me.cgrand.net/2009/10/15/multidim-arrays/
;; --------------------------------------------------------------------
    
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
         
;; --------------------------------------------------------------------
;; Main code
;; --------------------------------------------------------------------

(defn generate_table [n]
  (make-array Double/TYPE (inc n) (inc n)))
  
(defn flush_table [table]
  (do
    (dotimes [i (count table)]
      (dotimes [j (count table)]
        (deep-aset doubles table i j (if (and (<= i j) (<= j (+ i min_loop_size)) (>= i 1)) 1 0))))
    table))

(defn match_pairs_transient [rna_structure i j]
  (let [padded_structure (pad_string rna_structure)]
    (reduce 
      (fn [hash_of_pairs index]
        (do
          (cond
            (= \( (nth padded_structure index))
              (conj! (hash_of_pairs :stack) index)
            (= \) (nth padded_structure index))
              (if (-> (hash_of_pairs :stack) count zero? not)
                (do
                  (conj! (hash_of_pairs :set_of_pairs) [(nth (hash_of_pairs :stack) (-> (hash_of_pairs :stack) count dec)) index])
                  (pop! (hash_of_pairs :stack)))))
          hash_of_pairs))
      { :stack (transient []) :set_of_pairs (transient #{}) }
      (inc_range i j))))
      
(defn match_pairs [rna_structure i j]
  (if (> i j)
    #{}
    (-> (match_pairs_transient rna_structure i j) :set_of_pairs persistent!)))
    
(defn can_pair [rna_sequence i j]
  (let [padded_sequence (.toLowerCase (pad_string rna_sequence))]
    (some #{(nth padded_sequence j)} (pairings (nth padded_sequence i)))))

(defn end_paired [rna_structure i j]
  (some #{j} (map last (match_pairs rna_structure i j))))
      
(defn pair_distance [rna_structure i k j]
  (let [
    reference_structure (match_pairs rna_structure i j)
    comparitive_pairings (union (match_pairs rna_structure i (dec k)) (match_pairs rna_structure (inc k) (dec j)) #{[k j]})
    ] (count (symmetric_difference reference_structure comparitive_pairings))))
    
(defn partition_contribution [rna_structure table x_value i k j]
  (deep-aset doubles table i j
    (+
      (deep-aget doubles table i j)
      (*
        (Math/pow x_value (pair_distance rna_structure i k j))
        (if (= k i)
          (deep-aget doubles table (inc k) (dec j))
          (* 
            (deep-aget doubles table i (dec k)) 
            (deep-aget doubles table (inc k) (dec j))))))))

(defn solve_recurrences 
  ([rna_sequence rna_structure table x_value]
    (let [length (dec (count table))]
      (for [distance (inc_range (inc min_loop_size) (dec length)) i (inc_range 1 (- length distance))]
        (let [j (+ i distance)]
          (do
            (deep-aset doubles table i j 
              (* 
                (deep-aget doubles table i (dec j)) 
                (Math/pow x_value (if (end_paired rna_structure i j) 1 0))))
            (for [k (range i (- j min_loop_size)) :when (can_pair rna_sequence k j)]
              (partition_contribution rna_structure table x_value i k j))))))))
            
(def sample_sequence  "gggggccccc")
(def sample_structure "..........")
(def sample_table (flush_table (generate_table (count sample_sequence))))
(solve_recurrences sample_sequence sample_structure sample_table 1)
            
; (time (dotimes [_ 1000] (pair_distance ".(())(((())((.()..(.(.(.()().(.)()(.(.(.)((.().)..(((..)(.(.).(...(..(.()...)().)(().())).)..((.)..)((()).((.())).)().)(.).()).(.()))).).).((.().(.))(()()())(..)((.)).)(((((.(..()(.().(((().).)())().))))..(.)(.((())()()..)...)((..(.)).)).)((((..(()..(((.()..)()..).)())))(()))((.(.()()))((...(....)(().)((...(((.)()))..()))((...())).)))()().(.().()()(((..)....(((.()((..(.(..(.()().)).(..)))...)))..((...)(((()..((.)(.))....(.))((..)..(.((.)((()(((.)().)(.))(())..(.))()().(.())..)(..(...)(()(().)))..()..)))(()(().(((((()(()))(()))))...).)()(.)..))().())())()))(.)(.)(.((.(()(.(.(((().()...()()(.....(.((..)..)).).....))()()(.(.().(((((.()...).)...((.(((())()...))))(.()(()))))(..()))...().)(((.).))...(.)()()((.(.).(((.)).)(((..)..)((.()(.)())))...)))).)().(()))))((..)()((.(.().)(())).())))(())).)).((((()))).()()())))(.(.((..)...(.)))).))(..)..).()((.()(())()(((())(..).()(....).)(()()()))())(.)))....)()().).)))(.(..()...()((.)(.()()...(().()..)().(()).(.)..()().)).)()..(.(.)(.)).().)).).)))))." 1 500 1000)))
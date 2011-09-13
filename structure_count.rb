class Rnabor
  MIN_LOOP_SIZE = 3
  BASE_PAIRINGS = {
    "a" => %w[u],
    "u" => %w[a g],
    "g" => %w[c u],
    "c" => %w[g]
  }
  
  attr_reader :sequence, :structure
  
  def initialize(sequence, structure)
    @sequence, @structure = sequence, structure
  end
  
  def match_pairs(structure = structure)
    # So indicies go from 1..n
    structure[0, 0] = " "
    
    # Annotated structure hash from 0..length
    (@match_pairs ||= {})[structure] ||= structure.split(//).each_with_index.inject({}) do |hash, (symbol, index)|
      hash.tap do      
        case symbol
        when "(" then hash[index] = nil
        when ")" then hash[hash.select { |_, value| value.nil? }.keys.max] = index
        end
      end
    end
  end
  
  def can_pair?(i, j)
    # Assumes i, j are 1 indexed
    BASE_PAIRINGS[sequence[i - 1]].include?(sequence[j - 1])
  end
  
  def paired?(i, j)
    match_pairs[i] == j
  end
  
  def delta_neighbors(total_distance)
    table = generate_table(total_distance)
    
    (1..total_distance).each do |distance|
      ((MIN_LOOP_SIZE + 1)..(sequence.length - 1)).each do |base_pair_distance|
        (1..(sequence.length - base_pair_distance)).each do |i|
          j = i + base_pair_distance
          
          table[distance][i][j] = table[distance - (paired?(i, j) ? 1 : 0)][i][j - 1]
          
          (i..(j - MIN_LOOP_SIZE - 1)).each do |k|
            if can_pair?(k, j)
              # Saving these in variables because the structure is 0 indexed, so some conversions need to happen
              sectional_pairs          = match_pairs(structure[((i - 1)..(j - 1))]).values.reject(&:nil?).count
              sectional_pairs_before_k = match_pairs(structure[((i - 1)..(k - 2))]).values.reject(&:nil?).count
              sectional_pairs_after_k  = match_pairs(structure[(k..(j - 2))]).values.reject(&:nil?).count
              pairs_across_k           = sectional_pairs - sectional_pairs_before_k - sectional_pairs_after_k + (paired?(k, j) ? -1 : 1)
              
              # table[distance][i][j] += ...
            end
          end
        end
      end
    end
  end
  
  def generate_table(total_distance)
    (0..total_distance).map do |distance|
      # Future optimization, since i <= j, we only need a triangular of the matrix, so can just modify accessors for j to reduce space.
      (0..sequence.length).map do |i|
        (0..sequence.length).map do |j|
          if distance == 0 && i <= j
            1
          elsif distance > 0 && i <= j && j <= i + MIN_LOOP_SIZE
            0
          end
        end
      end
    end
  end
end

# rna = Rnabor.new("auacgccguaguau", "..(((...).))..")
# rna.delta_neighbors(2)
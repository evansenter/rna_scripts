require "../lagrange/lagrange.rb"

module Rnabor
  class Nussinov
    AVOGADROS_NUMBER   = 6.0221415e23
    BOLTZMANN_CONSTANT = 1.3806503e-23
    BOLTZMANN_FACTOR   = BOLTZMANN_CONSTANT * AVOGADROS_NUMBER
    TEMPERATURE        = 37
    BASE_PAIR_ENERGY   = Math::E ** (-1 / (BOLTZMANN_FACTOR * TEMPERATURE))
    MIN_LOOP_SIZE      = 3
    BASE_PAIRINGS      = {
      "a" => %w[u],
      "u" => %w[a g],
      "g" => %w[c u],
      "c" => %w[g]
    }

    attr_reader :length, :sequence, :structure, :table

    def initialize(sequence, structure)
      @length    = sequence.length
      @sequence  = (" " + sequence.downcase).freeze
      @structure = (" " + structure).freeze
      @table     = generate_table
    end

    def match_pairs(structure_to_match = structure)
      structure_to_match = " " + structure_to_match unless structure_to_match[0] == " "
      
      if structure_to_match.length > 2
        if structure_to_match == structure
          @paired_structure ||= get_pairings(structure)
        else
          get_pairings(structure_to_match)
        end
      else
        {}
      end
    end
    
    def get_pairings(structure)
      structure.split(//).each_with_index.inject({}) do |hash, (symbol, index)|
        hash.tap do      
          case symbol
          when "(" then hash[index] = nil
          when ")" then hash[hash.select { |_, value| value.nil? }.keys.max] = index
          end
        end
      end
    end

    def count_pairs(pair_hash)
      pair_hash.reject { |key, value| key.nil? || value.nil? }.count
    end

    def can_pair?(i, j)
      BASE_PAIRINGS[sequence[i]].include?(sequence[j])
    end

    def paired?(i, j)
      match_pairs[i] == j
    end
    
    def partition_function
      data = (0..length).map { |i| 1.0 / (i + 1) }.map do |x_value|
        [x_value, solve_recurrences(x_value)]
      end
      
      Lagrange.new(*data).coefficients
    end

    def solve_recurrences(x_value)
      flush_table
      
      ((MIN_LOOP_SIZE + 1)..(length - 1)).each do |base_pair_distance|
        (1..(length - base_pair_distance)).each do |i|
          j = i + base_pair_distance
          
          j_unpaired_contribution = table_at(i, j - 1) * x_value ** (paired?(i, j) ? 1 : 0)
          j_paired_contribution   = (i...j).select { |k| can_pair?(k, j) }.inject(0.0) do |sum, k|
            i_j_pairs                  = count_pairs(match_pairs(structure[i..j]))
            upstream_partition_pairs   = count_pairs(match_pairs(structure[i..(k - 1)]))
            downstream_partition_pairs = count_pairs(match_pairs(structure[(k + 1)..(j - 1)]))
            base_pair_difference       = i_j_pairs - upstream_partition_pairs - downstream_partition_pairs + (paired?(k, j) ? -1 : 1)
            
            BASE_PAIR_ENERGY * table_at(i, k - 1) * table_at(k + 1, j - 1) * x_value ** base_pair_difference
          end
          
          table_at(i, j, j_unpaired_contribution + j_paired_contribution)
        end
      end
      
      table_at(1, length)
    end

    def generate_table
      (0...length).map { Array.new(length) }
    end
    
    def flush_table
      (1..length).each do |i|
        (1..length).each do |j|
          table_at(i, j, j >= i + MIN_LOOP_SIZE ? 1.0 : 0.0)
        end
      end
    end
    
    def table_at(i, j, value = nil)
      value.nil? ? table[i - 1][j - 1] : table[i - 1][j - 1] = value
    end
  end
end

rna = Rnabor::Nussinov.new("acgccguaguacgccguagu", "(((...).))(((...).))")
rna.partition_function
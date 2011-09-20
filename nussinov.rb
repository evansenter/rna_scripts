require "set"
require "awesome_print"
require "../lagrange/lagrange.rb"

module Rnabor
  class Nussinov
    TEMPERATURE        = 37
    BOLTZMANN_CONSTANT = 0.0019872370936902486 # kcal / mol / K
    BASE_PAIR_ENERGY   = 0 # Math::E ** (1 / (BOLTZMANN_CONSTANT * (TEMPERATURE + 273.15)))
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
    
    def partition_function
      data = (0..length).map { |i| i }.map do |x_value|
        [x_value, solve_recurrences(x_value)]
      end
      
      # data = (0..length).map { |i| 1.0 / (i + 1) }.map do |x_value|
      #   [x_value, solve_recurrences(x_value)]
      # end
      
      Lagrange.new(*data).coefficients
    end

    def solve_recurrences(x_value)
      flush_table
      
      ((MIN_LOOP_SIZE + 1)..(length - 1)).each do |base_pair_distance|
        (1..(length - base_pair_distance)).each do |i|
          j = i + base_pair_distance
          
          j_unpaired_contribution = table_at(i, j - 1) * x_value ** (end_base_paired?(i, j) ? 1 : 0)
          j_paired_contribution   = (i..(j - MIN_LOOP_SIZE - 1)).select { |k| can_pair?(k, j) }.inject(0.0) do |sum, k|
            sum + BASE_PAIR_ENERGY * (k - 1 < i ? 1.0 : table_at(i, k - 1)) * table_at(k + 1, j - 1) * x_value ** pair_distance(i, k, j)
          end
          
          table_at(i, j, j_unpaired_contribution + j_paired_contribution)
        end
      end
      
      table_at(1, length)
    end
    
    def pair_distance(i, k, j)
      reference_structure  = closed_pairs(match_pairs(structure[i..j])).map { |from, to| [from + i - 1, to + i - 1] }.to_set
      upstream             = closed_pairs(match_pairs(k - 1 < i ? "" : structure[i..(k - 1)])).map { |from, to| [from + i - 1, to + i - 1] }
      downstream           = closed_pairs(match_pairs(structure[(k + 1)..(j - 1)])).map { |from, to| [from + k, to + k] }
      comparitive_pairings = (upstream + downstream + [[k, j]]).to_set
      
      ((reference_structure - comparitive_pairings) + (comparitive_pairings - reference_structure)).size
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
    
    def end_base_paired?(i, j)
      closed_pairs(match_pairs(structure[i..j])).values.compact.include?(j - i + 1)
    end

    def closed_pairs(pair_hash)
      pair_hash.reject { |key, value| key.nil? || value.nil? }
    end

    def count_pairs(pair_hash)
      closed_pairs(pair_hash).count
    end    

    def can_pair?(i, j)
      BASE_PAIRINGS[sequence[i]].include?(sequence[j])
    end

    def paired?(i, j)
      match_pairs[i] == j
    end

    def generate_table
      (0...length).map { Array.new(length) }
    end
    
    def flush_table
      (1..length).each do |i|
        (1..length).each do |j|
          table_at(i, j, i <= j && j <= i + MIN_LOOP_SIZE ? 1.0 : nil)
        end
      end
    end
    
    def table_at(i, j, value = nil)
      value.nil? ? table[i - 1][j - 1] : table[i - 1][j - 1] = value
    end
  end
end

rna = Rnabor::Nussinov.new(
  "ggggcccc", 
  "(.(...))"
)
rna.partition_function
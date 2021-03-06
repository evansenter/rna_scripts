require "awesome_print"
require "set"
require "complex"
require "benchmark"
require "../lagrange/lagrange.rb"

module Rnabor
  class Nussinov
    TEMPERATURE        = 37
    BOLTZMANN_CONSTANT = 0.0019872370936902486 # kcal / mol / K
    BASE_PAIR_ENERGY   = Math::E ** (1 / (BOLTZMANN_CONSTANT * (TEMPERATURE + 273.15)))
    MIN_LOOP_SIZE      = 3
    BASE_PAIRINGS      = {
      "a" => %w[u],
      "u" => %w[a g],
      "g" => %w[c u],
      "c" => %w[g]
    }

    attr_reader :length, :sequence, :structure, :table, :values

    def initialize(options = {})
      options         = { values: :roots_of_unity }.merge(options)
      @length         = options[:sequence].length
      @sequence       = (" " + options[:sequence].downcase).freeze
      @structure      = (" " + (structure || "." * length)).freeze
      @values         = options[:values]
      @table          = generate_table
    end
    
    def structure_count
      @unscaled_solutions = generate_x_values.map do |x_value|
        [x_value, solve_recurrences(x_value, 1)]
      end
      
      puts
      
      Lagrange.new(*@unscaled_solutions).coefficients
    end
    
    def partition_function
      @unscaled_solutions = generate_x_values.map do |x_value|
        [x_value, solve_recurrences(x_value, BASE_PAIR_ENERGY)]
      end
      
      puts
      
      coefficients = Lagrange.new(*@unscaled_solutions).coefficients
      
      (->(sum) { coefficients.map { |boltzmann_factor| boltzmann_factor / sum } })[coefficients.inject { |a, b| a + b }]
    end

    def solve_recurrences(x_value, energy)
      print "."
      
      flush_table
      
      ((MIN_LOOP_SIZE + 1)..(length - 1)).each do |distance|
        (1..(length - distance)).each do |i|
          j = i + distance
  
          table[i][j] = table[i][j - 1] * (x_value ** (end_base_paired?(i, j) ? 1 : 0))
          
          (i..(j - MIN_LOOP_SIZE - 1)).select { |k| can_pair?(k, j) }.each do |k|              
            base_pair_distance = pair_distance(i, k, j)
            
            if k == i
              table[i][j] += table[k + 1][j - 1] * energy * (x_value ** base_pair_distance)
            else
              table[i][j] += table[i][k - 1] * table[k + 1][j - 1] * energy * (x_value ** base_pair_distance)
            end
          end
        end
      end
        
      table[1][length]
    end
    
    def pair_distance(i, k, j)
      if (@memoized_pair_distance ||= {})[[i, k, j]]
        @memoized_pair_distance[[i, k, j]]
      else
        reference_structure  = match_pairs(structure[i..j]).map { |from, to| [from + i - 1, to + i - 1] }.to_set
        upstream             = match_pairs(k - 1 < i ? "" : structure[i..(k - 1)]).map { |from, to| [from + i - 1, to + i - 1] }
        downstream           = match_pairs(structure[(k + 1)..(j - 1)]).map { |from, to| [from + k, to + k] }
        comparitive_pairings = (upstream + downstream + [[k, j]]).to_set

        @memoized_pair_distance[[i, k, j]] = ((reference_structure - comparitive_pairings) + (comparitive_pairings - reference_structure)).size
      end
    end

    def match_all_pairs
      @memoized_match_all_pairs ||= get_pairings(structure)
    end

    def match_pairs(structure_to_match)
      structure_to_match = " " + structure_to_match unless structure_to_match[0] == " "
      
      if structure_to_match.length >= MIN_LOOP_SIZE + 3 # Padded space + minimum loop size + base pair on either side of loop
        get_pairings(structure_to_match)
      else
        {}
      end
    end
    
    def get_pairings(structure)
      stack = []
      
      structure.each_char.each_with_index.inject({}) do |hash, (symbol, index)|
        hash.tap do      
          case symbol
          when "(" then stack.push(index)
          when ")" then hash[stack.pop] = index unless stack.empty?
          end
        end
      end
    end
    
    def end_base_paired?(i, j)
      match_pairs(structure[i..j]).values.compact.include?(j - i + 1)
    end

    def can_pair?(i, j)
      BASE_PAIRINGS[sequence[i]].include?(sequence[j])
    end

    def paired?(i, j)
      match_all_pairs[i] == j
    end

    def generate_table
      (0..length).map { Array.new(length + 1) }
    end
    
    def flush_table
      (1..length).each do |i|
        (i..length).each do |j|
          table[i][j] = 1.0 if j <= i + MIN_LOOP_SIZE
        end
      end
    end
    
    def generate_x_values
      values == :roots_of_unity ? self.class.roots_of_unity(length) : values
    end
    
    def self.roots_of_unity(length, scaling = 1)
      (0..length).map do |i|
        Complex(scaling * Math.cos(2 * Math::PI * i / (length + 1)), scaling * Math.sin(2 * Math::PI * i / (length + 1)))
      end
    end
  end
end

# Benchmark.measure { (rna = Rnabor::Nussinov.new(sequence: ?g * 20 + ?c * 20)).partition_function }.real
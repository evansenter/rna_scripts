require "awesome_print"
require "set"
require "complex"
require "benchmark"
require "narray"

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

    attr_reader :length, :sequence, :structure, :table, :boltzmann_factor_table, :values

    def initialize(options = {})
      options         = { values: :roots_of_unity }.merge(options)
      @length         = options[:sequence].length
      @sequence       = (" " + options[:sequence].downcase).freeze
      @structure      = (" " + (structure || "." * length)).freeze
      @values         = options[:values]
      @table          = generate_table
    end
    
    def structure_count
      solve_with_energy(1)
    end
    
    def partition_function
      solve_with_energy(BASE_PAIR_ENERGY)
    end
    
    def solve_with_energy(energy)
      @unscaled_solutions = generate_x_values.map do |x_value|
        [x_value, solve_recurrences(x_value, energy)]
      end
      
      puts
      
      matrix_a = NMatrix[*@unscaled_solutions.map(&:first).map { |x| (0..length).map { |power| x ** power } }]
      vector_b = NVector[*@unscaled_solutions.map(&:last)]
      (vector_b / matrix_a).to_a
    end

    def solve_recurrences(x_value, energy)
      print "."
      
      flush_table
      
      ((MIN_LOOP_SIZE + 1)..(length - 1)).each do |distance|
        (1..(length - distance)).each do |i|
          j = i + distance
  
          table[i][j] = (table[i][j - 1] * (x_value ** end_base_paired?(i, j)))
          
          (i..(j - MIN_LOOP_SIZE - 1)).select { |k| can_pair?(k, j) }.each do |k|              
            base_pair_distance = 
            
            if k == i
              table[i][j] += table[k + 1][j - 1] * energy * x_value ** (base_pairs_matrix[i][j] - base_pairs_matrix[k + 1][j - 1] + paired?(k, j))
            else
              table[i][j] += table[i][k - 1] * table[k + 1][j - 1] * energy * (x_value ** (base_pairs_matrix[i][j] - base_pairs_matrix[i][k - 1] - base_pairs_matrix[k + 1][j - 1] + paired?(k, j)))
            end
          end
        end
      end
        
      table[1][length]
    end
    
    def get_pairings(structure)
      if instance_variable_defined?(:@base_pairings)
        @base_pairings
      else
        stack = []
      
        @base_pairings = structure.each_char.each_with_index.inject(Array.new(structure.length, -1)) do |array, (symbol, index)|
          array.tap do      
            case symbol
            when "(" then stack.push(index)
            when ")" then 
              if stack.empty?
                raise "Too many ')' in '#{structure}'"
              else
                stack.pop.tap do |opening|
                  array[opening] = index
                  array[index]   = opening
                end
              end
            end
          end
        end.tap do
          raise "Too many '(' in '#{structure}'" unless stack.empty?
        end
      end
    end
    
    def number_of_base_pairs(i, j)
      base_pairings = get_pairings(structure)
  
      (i..j).inject(0) do |count, index|
        count + (i < base_pairings[i] && j >= base_pairings[i] ? 1 : 0)
      end
    end

    def base_pairs_matrix
      @base_pairs_matrix ||= (0..length).map do |i|
        (0..length).map do |j|
          unless i.zero? || j.zero?
            number_of_base_pairs(i, j)
          end
        end
      end
    end
    
    def end_base_paired?(i, j)
      get_pairings(structure)[j] >= i ? 1 : 0
    end

    def can_pair?(i, j)
      BASE_PAIRINGS[sequence[i]].include?(sequence[j])
    end

    def paired?(i, j)
      get_pairings(structure)[i] == j ? -1 : 1
    end

    def generate_table
      (0..length).map { Array.new(length + 1, 1.0) }
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
    
    def self.roots_of_unity(length)
      (0..length).map do |i|
        Complex(Math.cos(2 * Math::PI * i / (length + 1)), Math.sin(2 * Math::PI * i / (length + 1)))
      end
    end
  end
end

if ARGV.empty?
  puts "Call: ruby ./nussinov_with_matricies.rb [SEQUENCE]"
else
  puts Rnabor::Nussinov.new(sequence: ARGV.first).structure_count
end
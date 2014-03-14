#!/hgsc_software/ruby/latest/bin/ruby
# -*- coding: utf-8 -*-

#name: SVachra_v1.9.rb
#author: Oliver A. Hampton
#date: Feb 25, 2014
#description: SVachra – Structural Variation Assesment of CHRomosomal Aberrations
#             A structural variation breakpoint caller that uses discordant mate pair reads consisting of both
#             inward and outward facing read types; such as the data delivered by Illumina mate pair and Nextera
#             sequencing libraries.  The SVachra breakpoint calling tool utilizes the characteristics of both
#             the inward and outward facing reads to call putative aberrant joints.  
#
#             NOTE: While inward and outward facing read types inform the assessment of structural variants, 
#                   reported SVs are orientated to the inward facing read orientations allowing easy comparison
#                   and integration of SV call with more common paired end sequencing data. 
#
#             SVachra is a play on the word “chakra” and refers to the ability to detect Structural Variation by
#             simultaneously evaluating inward and outward facing read data.
#
#program output:  Base_Name equals the bamfile name given by the user
#                 1. Base_Name.hist.txt         - Fragment_Lengths(bins) and Read_Pair_Counts to plot distributions of seq. library fragment sizes
#                 2. Base_Name.svp              - Main Output: listing of all structural variation annotations [types: INS,DEL,INV,ITX,CTX]
#                 3. Base_Name.bed              - Bed File of all structural variation annotations
#                 4. Base_Name.circos.link.txt  - Circos Link input file of all structural variation annotations
#

require 'getoptlong'
require 'fileutils'
require 'rubygems'
require 'mathn'

#CONSTANTS
#############################################################

##CHANGE LOCATION OF SAMTOOLS TO YOUR LOCATION
SAMTOOLS = "/hgsc_software/samtools/samtools-0.1.18/samtools"

##CONSTANTS
WINDOW = 100
INFINITY = 1.0/0
#K-MEANS CLUSTERING PARAMETERS
Knum = 3  
Sigma = 0
Delta = 0.001
#CLUSTER GROWTH PARAMETERS
SPAN = 2.5

#METHODS
#############################################################

def processArguments()
  optsArray =  [  ['--BAMFile',  '-f', GetoptLong::REQUIRED_ARGUMENT],
                  ['--ScreenOutBedFile', '-b', GetoptLong::OPTIONAL_ARGUMENT],
                  ['--min_cluster_count', '-c', GetoptLong::OPTIONAL_ARGUMENT],
                  ['--min_mapping_quality', '-m', GetoptLong::OPTIONAL_ARGUMENT],
                  ['--SVname', '-n', GetoptLong::OPTIONAL_ARGUMENT],
                  ['--SV_QC_Filtering', '-s', GetoptLong::NO_ARGUMENT],
                  ['--help', '-h', GetoptLong::NO_ARGUMENT]  ]

  progOpts = GetoptLong.new(*optsArray)
  missingOpts = Array.new
  optsArray.each{ |aa|
    if( aa[2] == 1 )
      missingOpts.push( aa[0] )
    end
  }
  
  optsHash = Hash.new{|hh,kk| hh[kk]=""}
  begin
    progOpts.each do |opt,arg|
      optsHash[opt] = arg
      missingOpts.delete( opt )
      end
    end

  if(missingOpts.length != 0 || optsHash.key?("--help") || optsHash.empty?)
    puts "SVachra_v1.9.rb -f reads.bam"
    puts "optional arguments:   -b screen out .bed file (optional)"
    puts "                      -c minimun number of non-overlapping read pairs to call breakpoint (DEFAULT=2)"
    puts "                      -m minimum mapping quaility threshold (DEFAULT=0)"
    puts "                      -n SV annoation name - appended with incrementing integer (DEFAULT=SV)"
    puts "                      -s SV cluster quality control filter based on proximity overlaps (optional)" 
    puts "help                  --help          [-h]"
  end
  return optsHash
end

##################################################################

module Enumerable
  
  def sum
    self.inject(0){|accum, i| accum + i }
  end
  
  def mean
    self.sum/self.length.to_f
  end
  
  def sample_variance
    m = self.mean
    sum = self.inject(0){|accum, i| accum +(i-m)**2 }
    (1/self.length.to_f*sum)
  end
  
  def standard_deviation
    return Math.sqrt(self.sample_variance)
  end
  
end 

########################################################################################### 
#K-MEANS Cluster class, represents a centroid point along with its associated nearby points

class KMeans_Cluster
  attr_accessor :center, :points

  #Constructor with a starting centerpoint                                                                                      
  def initialize(center)
    @center = center
    @points = []
  end

  #Recenters the centroid point and removes all of the associated points
  def recenter!
    xa = 0
    old_center = @center
    #Sum up all points
    @points.each do |point|
      xa += point
    end
    #Average out data 
    xa /= points.length
    #Reset center and return distance moved
    @center = xa
    return (old_center - center).abs
  end
end

#kmeans algorithm
def kmeans(data, k, delta=Delta)
  clusters = []
  data_unique = data.uniq
  #Assign intial values for all clusters
  (1..k).each do |point|
    if(data_unique.length == 0)
      puts "ERROR: INSERT SIZE DATA NOT UNIQUE"
      exit 1
    end 
    #Ruby 1.8 version implementation
    #rand_point = data.choice
    #Ruby 1.9 version implementation
    #rand_point = data.sample
    #Ruby 1.X version implementation
    index = (data_unique.length * rand).to_i
    rand_point = data_unique[index]
    c = KMeans_Cluster.new(rand_point)
    clusters.push c
    data_unique.delete(rand_point)
  end
  # Loop                                                        
  while true
    # Assign points to clusters
    data.each do |point|
      min_dist = +INFINITY
      min_cluster = nil
      # Find the closest cluster
      clusters.each do |cluster|
        dist = (point - cluster.center).abs        
        if dist < min_dist
          min_dist = dist
          min_cluster = cluster
        end
      end
      # Add to closest cluster
      min_cluster.points.push point
    end
    # Check deltas
    max_delta = -INFINITY
    clusters.each do |cluster|
      dist_moved = cluster.recenter!
      # Get largest delta 
      if dist_moved > max_delta
        max_delta = dist_moved
      end
    end
    # Check exit condition
    if max_delta < delta
      return clusters
    end
    # Reset points for the next iteration
    clusters.each do |cluster|
      cluster.points = []
    end
  end
end

##################################################################
#READ CLUSTER CLASS

class Cluster
  attr_accessor :chrm1, :pos1_min, :pos1_max, :ori1, :chrm2, :pos2_min, :pos2_max, :ori2, :read_id, :count, :size, :fragment, :indel, :type, :merge, :inv_merge, :inv_hash, :qc
  def initialize(line)
    parseBamLine(line)
  end

  def parseBamLine(theLine)
    arrSplit = theLine.split(/\t/)
    @read_id = [ arrSplit[0] ]
    @fragment = [ arrSplit[8].to_i.abs ]
    @indel = []
    @chrm1 = arrSplit[2]
    @pos1_min = arrSplit[3].to_i
    @pos1_max = arrSplit[3].to_i
    if( (arrSplit[1].to_i & 16) == 0 )
      @ori1 = '+'
    else
      @ori1 = '-'
    end
    if(arrSplit[6] == '=')
      @chrm2 = arrSplit[2]
    else
      @chrm2 = arrSplit[6]
    end
    @pos2_min = arrSplit[7].to_i
    @pos2_max = arrSplit[7].to_i
    if( (arrSplit[1].to_i & 32) == 0 )
      @ori2 = '+'
    else
      @ori2 = '-'
    end
    @count = 1
    @qc = 1
    @size = 0
    @merge = 0
    @inv_merge = 0
    @inv_hash = { "chrm" => nil, "ori" => nil, "index" => nil }
    @type = { "INS" => 0, "DEL" => 0, "INV" => 0, "ITX" => 0, "CTX" => 0, "UNK" => 0 }
    
    if( @chrm1 != @chrm2 )
      @type["CTX"] += 1
    else
      if( @ori1 == @ori2 )
        @type["INV"] += 1
      else
        if( arrSplit[3].to_i < arrSplit[7].to_i && @ori1 == "+" )
          @type["ITX"] += 1
        elsif( arrSplit[3].to_i < arrSplit[7].to_i && @ori1 == "-" )
          if( arrSplit[8].to_i.abs < PARAM_OUTWARD_MIN )
            @type["INS"] += 1
            @indel.push( PARAM_OUTWARD_MIN - arrSplit[8].to_i.abs )
          elsif( arrSplit[8].to_i.abs > PARAM_OUTWARD_MAX )
            @type["DEL"] += 1
            @indel.push( arrSplit[8].to_i.abs - PARAM_OUTWARD_MAX )
          else
            @type["UNK"] += 1
          end
        elsif( arrSplit[7].to_i < arrSplit[3].to_i && @ori2 == "+" )
          @type["ITX"] += 1
        elsif( arrSplit[7].to_i < arrSplit[3].to_i && @ori2 == "-" )
          if( arrSplit[8].to_i.abs < PARAM_OUTWARD_MIN )
            @type["INS"] +=1
            @indel.push( PARAM_OUTWARD_MIN - arrSplit[8].to_i.abs )
          elsif( arrSplit[8].to_i.abs > PARAM_OUTWARD_MAX )
            @type["DEL"] +=1
            @indel.push( arrSplit[8].to_i.abs - PARAM_OUTWARD_MAX )
          else
            @type["UNK"] += 1
          end
        end
      end
    end

  end


  def getType
    return self.type.max_by{|kk,vv| vv}[0]
  end


  def getOri(theInt)
    if( theInt == 1 )
      if( self.ori1 == "+" )
        return "-"
      elsif( self.ori1 == "-" )
        return "+"
      end
    elsif( theInt == 2 )
      if( self.ori2 == "+" )
        return "-"
      elsif( self.ori2 == "-" )
        return "+"
      end
    else
      return "ORI_ERROR"
    end
  end


  def intersect?(theLine)
    flag = 0
    arrSplit = theLine.split(/\t/)

    this_read_id = arrSplit[0]
    this_chrm1 = arrSplit[2]
    this_pos1 = arrSplit[3].to_i
    if( (arrSplit[1].to_i & 16) == 0 )
      this_ori1 = '+'
    else
      this_ori1 = '-'
    end
    if(arrSplit[6] == '=')
      this_chrm2 = arrSplit[2]
    else
      this_chrm2 = arrSplit[6]
    end
    this_pos2 = arrSplit[7].to_i
    if( (arrSplit[1].to_i & 32) == 0 )
      this_ori2 = '+'
    else
      this_ori2 = '-'
    end

    #CASE 0: SAME READ_ID
    if( self.read_id.include?( this_read_id ) )
      flag = 1
    else

      self_mid1 = (self.pos1_min + self.pos1_max)/2
      self_mid2 = (self.pos2_min + self.pos2_max)/2

      #CASE 1: SAME CHROMOSOME                                                                                                                            
      if( self.chrm1 == self.chrm2 && this_chrm1 == this_chrm2 && self.chrm1 == this_chrm1 )
        #SUBCASE A: SAME CHROMOSOME REARRANGEMENT; SELF_1 == LINE_1 && SELF_2 == LINE_2                                                                   
        if( ((self_mid1 - this_pos1).abs < (self_mid1 - this_pos2).abs) && ((self_mid2 - this_pos2).abs < (self_mid2 - this_pos1).abs) )
          if( self.ori1 == this_ori1 && self.ori2 == this_ori2 )
            if( self_mid1 >= (this_pos1 - PARAM_OUTWARD_MAX) && self_mid1 <= (this_pos1 + PARAM_OUTWARD_MAX) && self_mid2 >= (this_pos2 - PARAM_OUTWARD_MAX) && self_mid2 <= (this_pos2 + PARAM_OUTWARD_MAX) )
              flag = 1
            end
          end

          #SUBCASE B: SAME CHROMOSE REARRANGEMENT: SELF_1 == LINE_2 && SELF_2 == LINE_1                                                                   
        elsif( ((self_mid1 - this_pos2).abs < (self_mid1 - this_pos1).abs) && ((self_mid2 - this_pos1).abs < (self_mid2 - this_pos2).abs) )
          if( self.ori1 == this_ori2 && self.ori2 == this_ori1 )
            if( self_mid1 >= (this_pos2 - PARAM_OUTWARD_MAX) && self_mid1 <= (this_pos2 + PARAM_OUTWARD_MAX) && self_mid2 >= (this_pos1 - PARAM_OUTWARD_MAX) && self_mid2 <= (this_pos1 + PARAM_OUTWARD_MAX) )
              flag = 1
            end
          end
        end

        #CASE 2: DIFFERENT CHROMOSOME REARRANGEMENT; SELF_1 == LINE_1 && SELF_2 == LINE_2                                                                 
      elsif( self.chrm1 != self.chrm2 && self.chrm1 == this_chrm1 && self.chrm2 == this_chrm2 )
        if( self.ori1 == this_ori1 && self.ori2 == this_ori2 )
          if( self_mid1 >= (this_pos1 - PARAM_OUTWARD_MAX) && self_mid1 <= (this_pos1 + PARAM_OUTWARD_MAX) && self_mid2 >= (this_pos2 - PARAM_OUTWARD_MAX) && self_mid2 <= (this_pos2 + PARAM_OUTWARD_MAX) )
            flag = 1
          end
        end

        #CASE 3: DIFFERENT CHROMOSOME REARRANGEMENT; SELF_1 == LINE_2 && SELF_2 == LINE_1                                                                 
      elsif( self.chrm1 != self.chrm2 && self.chrm1 == this_chrm2 && self.chrm2 == this_chrm1 )
        if( self.ori1 == this_ori2 && self.ori2 == this_ori1 )
          if( self_mid1 >= (this_pos2 - PARAM_OUTWARD_MAX) && self_mid1 <= (this_pos2 + PARAM_OUTWARD_MAX) && self_mid2 >= (this_pos1 - PARAM_OUTWARD_MAX) && self_mid2 <= (this_pos1 + PARAM_OUTWARD_MAX) )
            flag = 1
          end
        end
      end
    end

    if( flag == 1 )
      return true
    else
      return false
    end
  end


  def cluster_intersect?(cluster)
    flag = 0
    #CASE 1: SAME CHROMOSOMES
    if( self.chrm1 == self.chrm2 && cluster.chrm1 == cluster.chrm2 && self.chrm1 == cluster.chrm1 )
      #SUBCASE A: SAME ORIENTATIONS
      if( self.ori1 == self.ori2 && cluster.ori1 == cluster.ori2 && self.ori1 == cluster.ori1 )
        if( (self.pos1_min <= cluster.pos1_max && cluster.pos1_min <= self.pos1_max) && (self.pos2_min <= cluster.pos2_max && cluster.pos2_min <= self.pos2_max) )
          flag = 1
        elsif( (self.pos1_min <= cluster.pos2_max && cluster.pos2_min <= self.pos1_max) && (self.pos2_min <= cluster.pos1_max && cluster.pos1_min <= self.pos2_max) )
          flag = 1
        end
        #SUBCASE B: DIFF ORIENTATIONS
      elsif( self.ori1 == cluster.ori1 && self.ori2 == cluster.ori2 )
        if( (self.pos1_min <= cluster.pos1_max&& cluster.pos1_min <= self.pos1_max) && (self.pos2_min <= cluster.pos2_max && cluster.pos2_min <= self.pos2_max) )
          flag = 1
        end
      elsif( self.ori1 == cluster.ori2 && self.ori2 == cluster.ori1 )
        if( (self.pos1_min <= cluster.pos2_max && cluster.pos2_min <= self.pos1_max) && (self.pos2_min <= cluster.pos1_max && cluster.pos1_min <= self.pos2_max) )
          flag = 1
        end
      end
      #CASE 2: DIFF CHROMOSOMES
    elsif( self.chrm1 == cluster.chrm1 && self.chrm2 == cluster.chrm2 )
      if( self.ori1 == cluster.ori1 && self.ori2 == cluster.ori2 )
        if( (self.pos1_min <= cluster.pos1_max&& cluster.pos1_min <= self.pos1_max) && (self.pos2_min <= cluster.pos2_max && cluster.pos2_min <= self.pos2_max) )
          flag = 1
        end
      end
    elsif( self.chrm1 == cluster.chrm2 && self.chrm2 == cluster.chrm1 )
      if( self.ori1 == cluster.ori2 && self.ori2 == cluster.ori1 )
        if( (self.pos1_min <= cluster.pos2_max && cluster.pos2_min <= self.pos1_max) && (self.pos2_min <= cluster.pos1_max && cluster.pos1_min <= self.pos2_max) )
          flag = 1
        end
      end
    end
    
   if( flag == 1 )
      return true
    else
      return false
    end
  end


  def addBP(theLine)
    arrSplit = theLine.split(/\t/)

    this_read_id = arrSplit[0]
    this_chrm1 = arrSplit[2]
    this_pos1 = arrSplit[3].to_i
    if( (arrSplit[1].to_i & 16) == 0 )
      this_ori1 = '+'
    else
      this_ori1 = '-'
    end
    if(arrSplit[6] == '=')
      this_chrm2 = arrSplit[2]
    else
      this_chrm2 = arrSplit[6]
    end
    this_pos2 = arrSplit[7].to_i
    if( (arrSplit[1].to_i & 32) == 0 )
      this_ori2 = '+'
    else
      this_ori2 = '-'
    end

    #CASE 0: SAME READ_ID
    if( self.read_id.include?( this_read_id ) )
      return true
    else

      self_mid1 = (self.pos1_min + self.pos1_max)/2
      self_mid2 = (self.pos2_min + self.pos2_max)/2

      #CASE 1: SAME CHROMOSOME
      if( self.chrm1 == self.chrm2 && this_chrm1 == this_chrm2 && self.chrm1 == this_chrm1 )
        #SUBCASE A: SAME CHROMOSOME REARRANGEMENT; SELF_1 == LINE_1 && SELF_2 == LINE_2
        if( ((self_mid1 - this_pos1).abs < (self_mid1 - this_pos2).abs) && ((self_mid2 - this_pos2).abs < (self_mid2 - this_pos1).abs) )
          if( self.ori1 == this_ori1 && self.ori2 == this_ori2 )
            if( self_mid1 >= (this_pos1 - PARAM_OUTWARD_MAX) && self_mid1 <= (this_pos1 + PARAM_OUTWARD_MAX) && self_mid2 >= (this_pos2 - PARAM_OUTWARD_MAX) && self_mid2 <= (this_pos2 + PARAM_OUTWARD_MAX) )
              if( this_pos1 < self.pos1_min )
                new_pos1_min = this_pos1
              else
                new_pos1_min = self.pos1_min
              end
              if( this_pos1 > self.pos1_max )
                new_pos1_max = this_pos1
              else
                new_pos1_max = self.pos1_max
              end
              if( this_pos2 < self.pos2_min )
                new_pos2_min = this_pos2
              else
                new_pos2_min = self.pos2_max
              end
              if( this_pos2 > self.pos2_max )
                new_pos2_max = this_pos2
              else
                new_pos2_max = self.pos2_max
              end
              if( (new_pos1_max - new_pos1_min) <= PARAM_OUTWARD_MAX && (new_pos2_max - new_pos2_min) <= PARAM_OUTWARD_MAX )
                self.pos1_min = new_pos1_min
                self.pos1_max = new_pos1_max
                self.pos2_min = new_pos2_min
                self.pos2_max = new_pos2_max
                self.size = (new_pos1_max - new_pos1_min) + (new_pos2_max - new_pos2_min)
                self.read_id.push( this_read_id )
                self.fragment.push( arrSplit[8].to_i.abs ) 
                if( this_chrm1 != this_chrm2 )
                  self.type["CTX"] += 1
                else
                  if( this_ori1 == this_ori2 )
                    self.type["INV"] += 1
                  else
                    if( this_pos1 < this_pos2 && this_ori1 == "+" )
                      self.type["ITX"] += 1
                    elsif( this_pos1 < this_pos2 && this_ori1 == "-" )
                      if( arrSplit[8].to_i.abs < PARAM_OUTWARD_MIN )
                        self.type["INS"] += 1
                        self.indel.push( PARAM_OUTWARD_MIN - arrSplit[8].to_i.abs )
                      elsif( arrSplit[8].to_i.abs > PARAM_OUTWARD_MAX )
                        self.type["DEL"] += 1
                        self.indel.push( arrSplit[8].to_i.abs - PARAM_OUTWARD_MAX )
                      else
                        self.type["UNK"] += 1
                      end
                    elsif( this_pos2 < this_pos1 && this_ori2 == "+" )
                      self.type["ITX"] += 1
                    elsif( this_pos2 < this_pos1 && this_ori2 == "-" )
                      if( arrSplit[8].to_i.abs < PARAM_OUTWARD_MIN )
                        self.type["INS"] += 1
                        self.indel.push( PARAM_OUTWARD_MIN - arrSplit[8].to_i.abs )
                      elsif( arrSplit[8].to_i.abs > PARAM_OUTWARD_MAX )
                        self.type["DEL"] += 1
                        self.indel.push( arrSplit[8].to_i.abs - PARAM_OUTWARD_MAX )
                      else
                        self.type["UNK"] += 1
                      end
                    end
                  end
                end
                self.count += 1
                return true
              else
                return false
              end
            end
          end

          #SUBCASE B: SAME CHROMOSE REARRANGEMENT: SELF_1 == LINE_2 && SELF_2 == LINE_1
        elsif( ((self_mid1 - this_pos2).abs < (self_mid1 - this_pos1).abs) && ((self_mid2 - this_pos1).abs < (self_mid2 - this_pos2).abs) )
          if( self.ori1 == this_ori2 && self.ori2 == this_ori1 )
            if( self_mid1 >= (this_pos2 - PARAM_OUTWARD_MAX) && self_mid1 <= (this_pos2 + PARAM_OUTWARD_MAX) && self_mid2 >= (this_pos1 - PARAM_OUTWARD_MAX) && self_mid2<= (this_pos1 + PARAM_OUTWARD_MAX) )
              if( this_pos1 < self.pos2_min )
                new_pos2_min = this_pos1
              else
                new_pos2_min = self.pos2_min
              end
              if( this_pos1 > self.pos2_max )
                new_pos2_max = this_pos1
              else
                new_pos2_max = self.pos2_max
              end
              if( this_pos2 < self.pos1_min )
                new_pos1_min = this_pos2
              else
                new_pos1_min = self.pos1_max
              end
              if( this_pos2 > self.pos1_max )
                new_pos1_max = this_pos2
              else
                new_pos1_max = self.pos1_max
              end
              if( (new_pos1_max - new_pos1_min) <= PARAM_OUTWARD_MAX && (new_pos2_max - new_pos2_min) <= PARAM_OUTWARD_MAX )
                self.pos1_min = new_pos1_min
                self.pos1_max = new_pos1_max
                self.pos2_min = new_pos2_min
                self.pos2_max = new_pos2_max
                self.size = (new_pos1_max - new_pos1_min) + (new_pos2_max - new_pos2_min)
                self.read_id.push( this_read_id )
                self.fragment.push( arrSplit[8].to_i.abs )
                if( this_chrm1 != this_chrm2 )
                  self.type["CTX"] += 1
                else
                  if( this_ori1 == this_ori2 )
                    self.type["INV"] += 1
                  else
                    if( this_pos1 < this_pos2 && this_ori1 == "+" )
                      self.type["ITX"] += 1
                    elsif( this_pos1 < this_pos2 && this_ori1 == "-" )
                      if( arrSplit[8].to_i.abs < PARAM_OUTWARD_MIN )
                        self.type["INS"] += 1
                        self.indel.push( PARAM_OUTWARD_MIN - arrSplit[8].to_i.abs )
                      elsif( arrSplit[8].to_i.abs > PARAM_OUTWARD_MAX )
                        self.type["DEL"] += 1
                        self.indel.push( arrSplit[8].to_i.abs - PARAM_OUTWARD_MAX )
                      else
                        self.type["UNK"] += 1
                      end
                    elsif( this_pos2 < this_pos1 && this_ori2 == "+" )
                      self.type["ITX"] += 1
                    elsif( this_pos2 < this_pos1 && this_ori2 == "-" )
                      if( arrSplit[8].to_i.abs < PARAM_OUTWARD_MIN )
                        self.type["INS"] += 1
                        self.indel.push( PARAM_OUTWARD_MIN - arrSplit[8].to_i.abs )
                      elsif( arrSplit[8].to_i.abs > PARAM_OUTWARD_MAX )
                        self.type["DEL"] += 1
                        self.indel.push( arrSplit[8].to_i.abs - PARAM_OUTWARD_MAX )
                      else
                        self.type["UNK"] += 1
                      end
                    end
                  end
                end
                self.count += 1
                return true
              else
                return false
              end
            end
          end
        end

        #CASE 2: DIFFERENT CHROMOSOME REARRANGEMENT; SELF_1 == LINE_1 && SELF_2 == LINE_2
      elsif( self.chrm1 != self.chrm2 && self.chrm1 == this_chrm1 && self.chrm2 == this_chrm2 )
        if( self.ori1 == this_ori1 && self.ori2 == this_ori2 )
          if( self_mid1 >= (this_pos1 - PARAM_OUTWARD_MAX) && self_mid1 <= (this_pos1 + PARAM_OUTWARD_MAX) && self_mid2 >= (this_pos2 - PARAM_OUTWARD_MAX) && self_mid2 <= (this_pos2 + PARAM_OUTWARD_MAX) )
            if( this_pos1 < self.pos1_min )
              new_pos1_min = this_pos1
            else
              new_pos1_min = self.pos1_min
            end
            if( this_pos1 > self.pos1_max )
              new_pos1_max = this_pos1
            else
              new_pos1_max = self.pos1_max
            end
            if( this_pos2 < self.pos2_min )
              new_pos2_min = this_pos2
            else
              new_pos2_min = self.pos2_max
            end
            if( this_pos2 > self.pos2_max )
              new_pos2_max = this_pos2
            else
              new_pos2_max = self.pos2_max
            end
            if( (new_pos1_max - new_pos1_min) <= PARAM_OUTWARD_MAX && (new_pos2_max - new_pos2_min) <= PARAM_OUTWARD_MAX )
              self.pos1_min = new_pos1_min
              self.pos1_max = new_pos1_max
              self.pos2_min = new_pos2_min
              self.pos2_max = new_pos2_max
              self.size = (new_pos1_max - new_pos1_min) + (new_pos2_max - new_pos2_min)
              self.read_id.push( this_read_id )
              self.fragment.push( arrSplit[8].to_i.abs )
              if( this_chrm1 != this_chrm2 )
                self.type["CTX"] += 1
              else
                if( this_ori1 == this_ori2 )
                  self.type["INV"] += 1
                else
                  if( this_pos1 < this_pos2 && this_ori1 == "+" )
                    self.type["ITX"] += 1
                  elsif( this_pos1 < this_pos2 && this_ori1 == "-" )
                    if( arrSplit[8].to_i.abs < PARAM_OUTWARD_MIN )
                      self.type["INS"] += 1
                      self.indel.push( PARAM_OUTWARD_MIN - arrSplit[8].to_i.abs )
                    elsif( arrSplit[8].to_i.abs > PARAM_OUTWARD_MAX )
                      self.type["DEL"] += 1
                      self.indel.push( arrSplit[8].to_i.abs - PARAM_OUTWARD_MAX )
                    else
                      self.type["UNK"] += 1
                    end
                  elsif( this_pos2 < this_pos1 && this_ori2 == "+" )
                    self.type["ITX"] += 1
                  elsif( this_pos2 < this_pos1 && this_ori2 == "-" )
                    if( arrSplit[8].to_i.abs < PARAM_OUTWARD_MIN )
                      self.type["INS"] += 1
                      self.indel.push( PARAM_OUTWARD_MIN - arrSplit[8].to_i.abs )
                    elsif( arrSplit[8].to_i.abs > PARAM_OUTWARD_MAX )
                      self.type["DEL"] += 1
                      self.indel.push( arrSplit[8].to_i.abs - PARAM_OUTWARD_MAX )
                    else
                      self.type["UNK"] += 1
                    end
                  end
                end
              end
              self.count += 1
              return true
            else
              return false
            end
          end
        end

        #CASE 3: DIFFERENT CHROMOSOME REARRANGEMENT; SELF_1 == LINE_2 && SELF_2 == LINE_1
      elsif( self.chrm1 != self.chrm2 && self.chrm1 == this_chrm2 && self.chrm2 == this_chrm1 )
        if( self.ori1 == this_ori2 && self.ori2 == this_ori1 )
          if( self_mid1 >= (this_pos2 - PARAM_OUTWARD_MAX) && self_mid1 <= (this_pos2 + PARAM_OUTWARD_MAX) && self_mid2 >= (this_pos1 - PARAM_OUTWARD_MAX) && self_mid2 <= (this_pos1 + PARAM_OUTWARD_MAX) )
            if( this_pos1 < self.pos2_min )
              new_pos2_min = this_pos1
            else
              new_pos2_min = self.pos2_min
            end
            if( this_pos1 > self.pos2_max )
              new_pos2_max = this_pos1
            else
              new_pos2_max = self.pos2_max
            end
            if( this_pos2 < self.pos1_min )
              new_pos1_min = this_pos2
            else
              new_pos1_min = self.pos1_max
            end
            if( this_pos2 > self.pos1_max )
              new_pos1_max = this_pos2
            else
              new_pos1_max = self.pos1_max
            end
            if( (new_pos1_max - new_pos1_min) <= PARAM_OUTWARD_MAX && (new_pos2_max - new_pos2_min) <= PARAM_OUTWARD_MAX )
              self.pos1_min = new_pos1_min
              self.pos1_max = new_pos1_max
              self.pos2_min = new_pos2_min
              self.pos2_max = new_pos2_max
              self.size = (new_pos1_max - new_pos1_min) + (new_pos2_max - new_pos2_min)
              self.read_id.push( this_read_id )
              self.fragment.push( arrSplit[8].to_i.abs )
              if( this_chrm1 != this_chrm2 )
                self.type["CTX"] += 1
              else
                if( this_ori1 == this_ori2 )
                  self.type["INV"] += 1
                else
                  if( this_pos1 < this_pos2 && this_ori1 == "+" )
                    self.type["ITX"] += 1
                  elsif( this_pos1 < this_pos2 && this_ori1 == "-" )
                    if( arrSplit[8].to_i.abs < PARAM_OUTWARD_MIN )
                      self.type["INS"] += 1
                      self.indel.push( PARAM_OUTWARD_MIN - arrSplit[8].to_i.abs )
                    elsif( arrSplit[8].to_i.abs > PARAM_OUTWARD_MAX )
                      self.type["DEL"] += 1
                      self.indel.push( arrSplit[8].to_i.abs - PARAM_OUTWARD_MAX )
                    else
                      self.type["UNK"] += 1
                    end
                  elsif( this_pos2 < this_pos1 && this_ori2 == "+" )
                    self.type["ITX"] += 1
                  elsif( this_pos2 < this_pos1 && this_ori2 == "-" )
                    if( arrSplit[8].to_i.abs < PARAM_OUTWARD_MIN )
                      self.type["INS"] += 1
                      self.indel.push( PARAM_OUTWARD_MIN - arrSplit[8].to_i.abs )
                    elsif( arrSplit[8].to_i.abs > PARAM_OUTWARD_MAX )
                      self.type["DEL"] += 1
                      self.indel.push( arrSplit[8].to_i.abs - PARAM_OUTWARD_MAX )
                    else
                      self.type["UNK"] += 1
                    end
                  end
                end
              end
              self.count += 1
              return true
            else
              return false
            end
          end
        end
      end
    end
  end


  def inward_outward_intersect?(cluster)
    flag = 0

    self_size = (self.pos1_max - self.pos1_min) + (self.pos2_max - self.pos2_min)
    self_mid1 = (self.pos1_min + self.pos1_max)/2
    self_mid2 = (self.pos2_min + self.pos2_max)/2
    cluster_size = (cluster.pos1_max - cluster.pos1_min) + (cluster.pos2_max - cluster.pos2_min)
    cluster_mid1 = (cluster.pos1_min + cluster.pos1_max)/2
    cluster_mid2 = (cluster.pos2_min + cluster.pos2_max)/2
 
    if( self_size != cluster_size )
      #CASE 1: SAME CHROMOSOME
      if( self.chrm1 == self.chrm2 && cluster.chrm1 == cluster.chrm2 && self.chrm1 == cluster.chrm1 )
        #SUBCASE A: SAME CHROMOSOME REARRANGEMENT; SELF_1 == CLUSTER_1 && SELF_2 == CLUSTER_2
        if( ((self_mid1 - cluster_mid1).abs < (self_mid1 - cluster_mid2).abs) && ((self_mid2 - cluster_mid2).abs < (self_mid2 - cluster_mid2).abs) )
          if( self.ori1 != cluster.ori1 && self.ori2 != cluster.ori2 )
            if( self_mid1 >= (cluster_mid1 - PARAM_OUTWARD_MAX) && self_mid1 <= (cluster_mid1 + PARAM_OUTWARD_MAX) && self_mid2 >= (cluster_mid2 - PARAM_OUTWARD_MAX) && self_mid2 <= (cluster_mid2 + PARAM_OUTWARD_MAX) )
              flag = 1
            end
          end

        #SUBCASE B: SAME CHROMOSE REARRANGEMENT: SELF_1 == CLUSTER_2 && SELF_2 == CLUSTER_1
        elsif( ((self_mid1 - cluster_mid2).abs < (self_mid1 - cluster_mid1).abs) && ((self_mid2 - cluster_mid1).abs < (self_mid2 - cluster_mid2).abs) )
          if( self.ori1 != cluster.ori2 && self.ori2 != cluster.ori1 )
            if( self_mid1 >= (cluster_mid2 - PARAM_OUTWARD_MAX) && self_mid1 <= (cluster_mid2 + PARAM_OUTWARD_MAX) && self_mid2 >= (cluster_mid1 - PARAM_OUTWARD_MAX) && self_mid2<= (cluster_mid1 + PARAM_OUTWARD_MAX) )
              flag = 1
            end
          end
        end

        #CASE 2: DIFFERENT CHROMOSOME REARRANGEMENT; SELF_1 == CLUSTER_1 && SELF_2 == CLUSTER_2
      elsif( self.chrm1 != self.chrm2 && self.chrm1 == cluster.chrm1 && self.chrm2 == cluster.chrm2 )
        if( self.ori1 == cluster.ori1 && self.ori2 == cluster.ori2 )
          if( self_mid1 >= (cluster_mid1 - PARAM_OUTWARD_MAX) && self_mid1 <= (cluster_mid1 + PARAM_OUTWARD_MAX) && self_mid2 >= (cluster_mid2 - PARAM_OUTWARD_MAX) && self_mid2 <= (cluster_mid2 + PARAM_OUTWARD_MAX) )
            flag = 1
          end
        end

      #CASE 3: DIFFERENT CHROMOSOME REARRANGEMENT; SELF_1 == CLUSTER_2 && SELF_2 == CLUSTER_1
      elsif( self.chrm1 != self.chrm2 && self.chrm1 == cluster.chrm2 && self.chrm2 == cluster.chrm1 )
        if( self.ori1 != cluster.ori2 && self.ori2 == cluster.ori1 )
          if( self_mid1 >= (cluster_mid2 - PARAM_OUTWARD_MAX) && self_mid1 <= (cluster_mid2 + PARAM_OUTWARD_MAX) && self_mid2 >= (cluster_mid1 - PARAM_OUTWARD_MAX) && self_mid2 <= (cluster_mid1 + PARAM_OUTWARD_MAX) )
            flag = 1
          end
        end
      end
    end

    if( flag == 1 )
      return true
    else
      return false
    end
  end

  def inward_outward_merge(cluster)
    self_size = (self.pos1_max - self.pos1_min) + (self.pos2_max - self.pos2_min)
    self_mid1 = (self.pos1_min + self.pos1_max)/2
    self_mid2 = (self.pos2_min + self.pos2_max)/2
    cluster_size = (cluster.pos1_max - cluster.pos1_min) + (cluster.pos2_max - cluster.pos2_min)
    cluster_mid1 = (cluster.pos1_min + cluster.pos1_max)/2
    cluster_mid2 = (cluster.pos2_min + cluster.pos2_max)/2

    #CASE 1: SAME CHROMOSOME
    if( self.chrm1 == self.chrm2 && cluster.chrm1 == cluster.chrm2 && self.chrm1 == cluster.chrm1 )
      #SUBCASE A: SAME CHROMOSOME REARRANGEMENT; SELF_1 == CLUSTER_1 && SELF_2 == CLUSTER_2
      if( ((self_mid1 - cluster_mid1).abs < (self_mid1 - cluster_mid2).abs) && ((self_mid2 - cluster_mid2).abs < (self_mid2 - cluster_mid1).abs) )
        if( self.ori1 != cluster.ori1 && self.ori2 != cluster.ori2 )
          if( self_mid1 >= (cluster_mid1 - PARAM_OUTWARD_MAX) && self_mid1 <= (cluster_mid1 + PARAM_OUTWARD_MAX) && self_mid2 >= (cluster_mid2 - PARAM_OUTWARD_MAX) && self_mid2 <= (cluster_mid2 + PARAM_OUTWARD_MAX) )

            if( cluster.pos1_min < self.pos1_min )
              new_pos1_min = cluster.pos1_min
            else
              new_pos1_min = self.pos1_min
            end
            if( cluster.pos1_max > self.pos1_max )
              new_pos1_max = cluster.pos1_max
            else
              new_pos1_max = self.pos1_max
            end
            if( cluster.pos2_min < self.pos2_min )
              new_pos2_min = cluster.pos2_min
            else
              new_pos2_min = self.pos2_max
            end
            if( cluster.pos2_max > self.pos2_max )
              new_pos2_max = cluster.pos2_max
            else
              new_pos2_max = self.pos2_max
            end

            if( ((new_pos1_max - new_pos1_min) + (new_pos2_max - new_pos2_min)) <= (PARAM_OUTWARD_MAX * SPAN).to_i )
              self.pos1_min = new_pos1_min
              self.pos1_max = new_pos1_max
              self.pos2_min = new_pos2_min
              self.pos2_max = new_pos2_max
              self.size = (new_pos1_max - new_pos2_min) + (new_pos2_max - new_pos2_min)
              self.read_id = self.read_id + cluster.read_id
              self.fragment = self.fragment + cluster.fragment
              self.count = self.count + cluster.count
              if( cluster_size > self_size )
                self.ori1 = cluster.ori1
                self.ori2 = cluster.ori2
                self.type = cluster.type
              end
              self.merge = 1
              return true
            else
              return false
            end
          end
        end

        #SUBCASE B: SAME CHROMOSE REARRANGEMENT: SELF_1 == CLUSTER_2 && SELF_2 == CLUSTER_1
      elsif( ((self_mid1 - cluster_mid2).abs < (self_mid1 - cluster_mid1).abs) && ((self_mid2 - cluster_mid1).abs < (self_mid2 - cluster_mid2).abs) )
        if( self.ori1 != cluster.ori2 && self.ori2 != cluster.ori1 )
          if( self_mid1 >= (cluster_mid2 - PARAM_OUTWARD_MAX) && self_mid1 <= (cluster_mid2 + PARAM_OUTWARD_MAX) && self_mid2 >= (cluster_mid1 - PARAM_OUTWARD_MAX) && self_mid2<= (cluster_mid1 + PARAM_OUTWARD_MAX) )

            if( cluster.pos1_min < self.pos2_min )
              new_pos2_min = cluster.pos1_min
            else
              new_pos2_min = self.pos2_min
            end
            if( cluster.pos1_max > self.pos2_max )
              new_pos2_max = cluster.pos1_max
            else
              new_pos2_max = self.pos2_max
            end
            if( cluster.pos2_min < self.pos1_min )
              new_pos1_min = cluster.pos2_min
            else
              new_pos1_min = self.pos1_max
            end
            if( cluster.pos2_max > self.pos1_max )
              new_pos1_max = cluster.pos2_max
            else
              new_pos1_max = self.pos1_max
            end

            if( ((new_pos1_max - new_pos1_min) + (new_pos2_max - new_pos2_min)) <= (PARAM_OUTWARD_MAX * SPAN).to_i )
              self.pos1_min = new_pos1_min
              self.pos1_max = new_pos1_max
              self.pos2_min = new_pos2_min
              self.pos2_max = new_pos2_max
              self.size = (new_pos1_max - new_pos2_min) + (new_pos2_max - new_pos2_min)
              self.read_id = self.read_id + cluster.read_id
              self.fragment = self.fragment + cluster.fragment
              self.count = self.count + cluster.count
              if( cluster_size > self_size )
                self.ori1 = cluster.ori2
                self.ori2 = cluster.ori1
                self.type = cluster.type
              end
              self.merge = 1
              return true
            else
              return false
            end
          end
        end
      end

      #CASE 2: DIFFERENT CHROMOSOME REARRANGEMENT; SELF_1 == CLUSTER_1 && SELF_2 == CLUSTER_2
    elsif( self.chrm1 != self.chrm2 && self.chrm1 == cluster.chrm1 && self.chrm2 == cluster.chrm2 )
      if( self.ori1 != cluster.ori1 && self.ori2 != cluster.ori2 )
        if( self_mid1 >= (cluster_mid1 - PARAM_OUTWARD_MAX) && self_mid1 <= (cluster_mid1 + PARAM_OUTWARD_MAX) && self_mid2 >= (cluster_mid2 - PARAM_OUTWARD_MAX) && self_mid2 <= (cluster_mid2 + PARAM_OUTWARD_MAX) )

          if( cluster.pos1_min < self.pos1_min )
            new_pos1_min = cluster.pos1_min
          else
            new_pos1_min = self.pos1_min
          end
          if( cluster.pos1_max > self.pos1_max )
            new_pos1_max = cluster.pos1_max
          else
            new_pos1_max = self.pos1_max
          end
          if( cluster.pos2_min < self.pos2_min )
            new_pos2_min = cluster.pos2_min
          else
            new_pos2_min = self.pos2_max
          end
          if( cluster.pos2_max > self.pos2_max )
            new_pos2_max = cluster.pos2_max
          else
            new_pos2_max = self.pos2_max
          end

          if( ((new_pos1_max - new_pos1_min) + (new_pos2_max - new_pos2_min)) <= (PARAM_OUTWARD_MAX * SPAN).to_i )
            self.pos1_min = new_pos1_min
            self.pos1_max = new_pos1_max
            self.pos2_min = new_pos2_min
            self.pos2_max = new_pos2_max
            self.size = (new_pos1_max - new_pos2_min) + (new_pos2_max - new_pos2_min)
            self.read_id = self.read_id + cluster.read_id
            self.fragment = self.fragment + cluster.fragment
            self.count = self.count + cluster.count
            if( cluster_size > self_size )
              self.ori1 = cluster.ori1
              self.ori2 = cluster.ori2
              self.type = cluster.type
            end
            self.merge = 1
            return true
          else
            return false
          end
        end
      end

    #CASE 3: DIFFERENT CHROMOSOME REARRANGEMENT; SELF_1 == CLUSTER_2 && SELF_2 == CLUSTER_1
    elsif( self.chrm1 != self.chrm2 && self.chrm1 == cluster.chrm2 && self.chrm2 == cluster.chrm1 )
      if( self.ori1 != cluster.ori2 && self.ori2 != cluster.ori1 )
        if( self_mid1 >= (cluster_mid2 - PARAM_OUTWARD_MAX) && self_mid1 <= (cluster_mid2 + PARAM_OUTWARD_MAX) && self_mid2 >= (cluster_mid1 - PARAM_OUTWARD_MAX) && self_mid2 <= (cluster_mid1 + PARAM_OUTWARD_MAX) )

          if( cluster.pos1_min < self.pos2_min )
            new_pos2_min = cluster.pos1_min
          else
            new_pos2_min = self.pos2_min
          end
          if( cluster.pos1_max > self.pos2_max )
            new_pos2_max = cluster.pos1_max
          else
            new_pos2_max = self.pos2_max
          end
          if( cluster.pos2_min < self.pos1_min )
            new_pos1_min = cluster.pos2_min
          else
            new_pos1_min = self.pos1_max
          end
          if( cluster.pos2_max > self.pos1_max )
            new_pos1_max = cluster.pos2_max
          else
            new_pos1_max = self.pos1_max
          end

          if( ((new_pos1_max - new_pos1_min) + (new_pos2_max - new_pos2_min)) <= (PARAM_OUTWARD_MAX * SPAN).to_i )
            self.pos1_min = new_pos1_min
            self.pos1_max = new_pos1_max
            self.pos2_min = new_pos2_min
            self.pos2_max = new_pos2_max
            self.size = (new_pos1_max - new_pos2_min) + (new_pos2_max - new_pos2_min)
            self.read_id = self.read_id + cluster.read_id
            self.fragment = self.fragment + cluster.fragment
            self.count = self.count + cluster.count
            if( cluster_size > self_size )
              self.ori1 = cluster.ori2
              self.ori2 = cluster.ori1
              self.type = cluster.type
            end
            self.merge = 1
            return true
          else
            return false
          end
        end
      end
    end
  end


  def inverse_intersect?(cluster,kk,oo,ii,jj)
    flag = 0

    if( self.chrm1 == self.chrm2 && cluster.chrm1 == cluster.chrm2 && self.chrm1 == cluster.chrm1 )
      if( self.ori1 == self.ori2 && cluster.ori1 == cluster.ori2 && self.ori1 != cluster.ori1 )
        if( ((self.pos1_min - PARAM_OUTWARD_MAX) <= (cluster.pos1_max + PARAM_OUTWARD_MAX) && (cluster.pos1_min - PARAM_OUTWARD_MAX) <= (self.pos1_max + PARAM_OUTWARD_MAX)) && ((self.pos2_min - PARAM_OUTWARD_MAX) <= (cluster.pos2_max + PARAM_OUTWARD_MAX) && (cluster.pos2_min - PARAM_OUTWARD_MAX) <= (self.pos2_max + PARAM_OUTWARD_MAX)) )
          flag = 1
        elsif( ((self.pos1_min - PARAM_OUTWARD_MAX) <= (cluster.pos2_max + PARAM_OUTWARD_MAX) && (cluster.pos2_min - PARAM_OUTWARD_MAX) <= (self.pos1_max + PARAM_OUTWARD_MAX)) && ((self.pos2_min - PARAM_OUTWARD_MAX) <= (cluster.pos1_max + PARAM_OUTWARD_MAX) && (cluster.pos1_min - PARAM_OUTWARD_MAX) <= (self.pos2_max + PARAM_OUTWARD_MAX)) )
          flag = 1
        end
      end
    end
    
    if( flag == 1 )
      self.inv_merge = 1
      self.inv_hash["chrm"] = kk
      self.inv_hash["ori"] = oo
      self.inv_hash["index"] = jj
      cluster.inv_merge = 1
      cluster.inv_hash["chrm"] = kk
      cluster.inv_hash["ori"] = oo
      cluster.inv_hash["index"] = ii
      return true
    else
      return false
    end
  end

end

##################################################################

optHash = processArguments()

screen_out_hash = Hash.new{|hh,kk| hh[kk]=Array.new}
if( optHash.key?("--ScreenOutBedFile") )
  screen_out_file = File.open( optHash["--ScreenOutBedFile"] )
  screen_out_file.each_line{ |line|
    line.chomp!
    arrSplit = line.split(/\t/)
    screen_out_hash[arrSplit[0]].push("#{arrSplit[1]}-#{arrSplit[2]}")
  }
end

if( optHash.key?("--min_mapping_quality") )
  MIN_MAPPING_QUALITY = optHash["--min_mapping_quality"].to_i
else
  MIN_MAPPING_QUALITY = 0
end

puts "STARTING: read in #{optHash["--BAMFile"]} ... accumulating insert sizes"

cmd = `#{SAMTOOLS} view -f1 -F1804 #{optHash["--BAMFile"]}`

counter = 0
total = cmd.count("\n")
puts "TOTAL BAM FILE ENTRIES: #{total}"

base_name = optHash["--BAMFile"]

hist_out = base_name + ".hist.txt"
hist_out_file = File.open(hist_out, "w") 

hist_hash = Hash::new(0)
unique_read_id = Hash::new(0)
seq_length = 0


cmd.each_line{ |line|
  counter += 1
  if( counter % 1000000 == 0 )
    print "#{counter/1000000}M ( #{((counter/total.to_f)*100).to_i}% ) ... "
  end 

  line.chomp!
  arrSplit = line.split(/\t/)
  if( arrSplit[9].length > seq_length )
    seq_length = arrSplit[9].length
  end
  if( arrSplit[4].to_i >= MIN_MAPPING_QUALITY )
    pass_screen = true
    screen_out_hash[arrSplit[2]].each{ |cc|
      arrCC = cc.split(/-/)
      if( (arrCC[0].to_i <= arrSplit[3].to_i && arrSplit[3].to_i <= arrCC[1].to_i) || (arrCC[0].to_i <= arrSplit[7].to_i && arrSplit[7].to_i <= arrCC[1].to_i) )
        pass_screen = false
      end
    }
    if( pass_screen == true )
      unique_read_id[arrSplit[0]] += 1
    end

    if( arrSplit[6] == '=' && pass_screen == true )
      if( unique_read_id[arrSplit[0]] == 2 )
        key = ( arrSplit[8].to_i.abs / WINDOW ).to_i
        hist_hash[key] += 1
      end
    end
  end
}      

hist_out_file.puts "Fragment_Lenghts(bins)\tRead_Pair_Counts"
hist_hash.keys.sort.each{|kk|
  hist_out_file.puts "#{kk.to_i*100}\t#{hist_hash[kk]}"
}
hist_out_file.close

counter = 0
print "done\n"

#############################################

puts "FINISHED: read in #{optHash["--BAMFile"]} file for insert sizes"

puts "STARTING: K-MEANS CLUSTERING of INSERT SIZE FREQUENCIES"

clusters = kmeans( hist_hash.values.uniq, Knum )

puts "FINISHED: K-MEANS CLUSTERING"

##DEBUG
#INSPECT KMEAN CLUSTERS
#clusters.each_index{ |ii|
#  puts "Uniq Cluster #{ii} Size:  #{clusters[ii].points.size}"
#  puts "Uniq Cluster #{ii} Min:   #{clusters[ii].points.min}"
#  puts "Uniq Cluster #{ii} Mean:  #{clusters[ii].points.mean}"
#  puts "Uniq Cluster #{ii} SD:    #{clusters[ii].points.standard_deviation}"
#  puts "Uniq Cluster #{ii} Max:   #{clusters[ii].points.max}"
#}                        
##

min_index = nil
min_mean = INFINITY
clusters.each_index{ |ii|
  if( clusters[ii].points.mean < min_mean )
    min_mean = clusters[ii].points.mean
    min_index = ii
  end
}

sigma = Sigma
background_noise = true

while( background_noise == true && sigma <= 3 )
  background_noise = false
  cutoff = ((clusters[min_index].points.mean)+(sigma*clusters[min_index].points.standard_deviation)).to_i
  ##SHOW BACKGROUND NOISE THRESHOLD                                                                          
  #puts "Sigma(#{sigma}): BACKGROUND CUTOFF (frequency): #{cutoff}"

  max_hist_index = hist_hash.key(hist_hash.values.max)
  if( hist_hash[max_hist_index] < cutoff )
    background_noise = true
    puts "BACKGROUND NOISE WITHIN DISTRIBUTION OF FRAGMENTS TOO HIGH. INCREASING SIGMA VALUE."
  end

  outward_min = nil
  outward_max = nil

  hist_index_arr = hist_hash.keys.sort

  outward_min_arr = hist_index_arr.slice(hist_index_arr.first..max_hist_index).reverse
  outward_min_arr.each_index{ |ii|
    if( hist_hash[outward_min_arr[ii]] < cutoff && outward_min.nil? )
      unless( outward_min_arr[ii-1].nil? )
        outward_min = outward_min_arr[ii-1]
      else
        outward_min = outward_min_arr[ii]
      end
    end
  }
  if(outward_min.nil?)
    outward_min = outward_min_arr.last
  end
  
  outward_max_arr = hist_index_arr.slice(max_hist_index..hist_index_arr.last)
  outward_max_arr.each_index{ |ii|
    if( hist_hash[outward_max_arr[ii]] < cutoff && outward_max.nil? )
      unless( outward_max_arr[ii-1].nil? )
        outward_max = outward_max_arr[ii-1]
      else
        outward_max = outward_max_arr[ii]
      end
    end 
  }
  if(outward_max.nil?)
    outward_max = outward_max_arr.last
  end
  
  second_hist_index_arr = hist_index_arr.dup
  hist_index_arr.each{ |ii|
    if( ii.to_i >= outward_min.to_i && ii.to_i <= outward_max.to_i )
      second_hist_index_arr.delete(ii)
    end
  }
  
  second_hist_hash = Hash.new
  second_hist_index_arr.each{ |ii|
    second_hist_hash[ii] = hist_hash[ii]
  }
  
  second_max_hist_index = second_hist_hash.key(second_hist_hash.values.max)
  if( second_hist_hash[second_max_hist_index] < cutoff )
    background_noise = true
    puts "BACKGROUND NOISE WITHIN DISTRIBUTION OF FRAGMENTS TOO HIGH. INCREASING SIGMA VALUE." 
  end

  inward_min = nil
  inward_max = nil

  second_hist_index_arr = second_hist_hash.keys.sort
  
  inward_min_arr = second_hist_index_arr.slice(second_hist_index_arr.first..second_max_hist_index).reverse
  inward_min_arr.each_index{ |ii|
    if( second_hist_hash[inward_min_arr[ii]] < cutoff && inward_min.nil? )
      unless( inward_min_arr[ii-1].nil? )
        inward_min = inward_min_arr[ii-1]
      else
        inward_min = inward_min_arr[ii]
      end
    end
  }
  if(inward_min.nil?)
    inward_min = inward_min_arr.last
  end
  
  inward_max_arr = second_hist_index_arr.slice(second_max_hist_index..second_hist_index_arr.last)
  inward_max_arr.each_index{ |ii|
    if( second_hist_hash[inward_max_arr[ii]] < cutoff && inward_max.nil? )
      unless( inward_max_arr[ii-1].nil? )
        inward_max = inward_max_arr[ii-1]
      else
        inward_max = inward_max_arr[ii]
      end
    end
  }
  if(inward_max.nil?)
    inward_max = inward_max_arr.last
  end

  sigma += 1
end

if( background_noise )
  puts "ERROR: POOR REPRESENTATION OF INSERT LIBRARY - BREAKPOINT CALLING TERMINATED"                                                                                                           
  exit 1     
end

unless( outward_min <= inward_max && inward_min <= outward_max )
  if( inward_max > outward_max )
    tmp_inward_min = inward_min
    tmp_inward_max = inward_max
    inward_min = outward_min
    inward_max = outward_max
    outward_min = tmp_inward_min
    outward_max = tmp_inward_max
  end
else
  puts "ERROR: INWARD and OUTWARD INSERT DISTRIBUTIONS NOT DISCRETE ENOUGH TO DECONVOLUTE -- UNABLE TO CALL STUCTURAL VARIATIONS"
  exit 1
end

clusters.clear
hist_index_arr.clear
hist_hash.clear
second_hist_hash.clear
second_hist_index_arr.clear

#INPUT PARAMETERS FOR INWARD AND OUTWARD FRAGMENT DISTRUBUTIONS
inward_min =  inward_min*WINDOW
inward_max = inward_max*WINDOW
outward_min = outward_min*WINDOW
outward_max = outward_max*WINDOW

unless( inward_min == 0 )
  inward_min = inward_min - WINDOW
end
inward_max = inward_max + WINDOW
unless( outward_min == 0 )
  outward_min = outward_min - WINDOW
end
outward_max = outward_max + WINDOW
outward_mid = (outward_min + outward_max)/2

puts "parameter - INWARD_MIN: #{inward_min}"
puts "parameter - INWARD_MAX: #{inward_max}"
puts "parameter - OUTWARD_MIN: #{outward_min}"
puts "parameter - OUTWARD_MAX: #{outward_max}"

##################################################################################

PARAM_INWARD_MIN = inward_min
PARAM_INWARD_MAX = inward_max
PARAM_OUTWARD_MIN = outward_min
PARAM_OUTWARD_MAX = outward_max

if( optHash.key?("--min_cluster_count") )
  MIN_CLUSTER_COUNT = optHash["--min_cluster_count"].to_i
else
  MIN_CLUSTER_COUNT = 2
end

clusterHash = Hash.new{|hh,kk| hh[kk]=Hash.new{|mm,nn| mm[nn]=Array.new}}

puts "STARTING: Reading in #{optHash["--BAMFile"]} BAM and Building Discordant Read Clusters..."

pass_counter = 0
cmd.each_line{ |line|
  counter += 1
  line.chomp!
  arrSplit = line.split(/\t/)
  if( arrSplit[4].to_i >= MIN_MAPPING_QUALITY )
    id = arrSplit[0]
    chrm1 = arrSplit[2]
    pos1 = arrSplit[3].to_i
    if(arrSplit[6] == '=')
      chrm2 = arrSplit[2]
    else
      chrm2 = arrSplit[6]
    end
    pos2 = arrSplit[7].to_i
    if( (arrSplit[1].to_i & 16) == 0 )
      ori1 = '+'
    else
      ori1 = '-'
    end
    if( (arrSplit[1].to_i & 32) == 0 )
      ori2 = '+'
    else
      ori2 = '-'
    end
    
    pass_screen = true
    ##SCREEN OUT BED INTERSECTIONS
    screen_out_hash[chrm1].each{ |cc|
      arrCC = cc.split(/-/)
      if( arrCC[0].to_i <= pos1 && pos1 <= arrCC[1].to_i )
        pass_screen = false
      end
    }
    screen_out_hash[chrm2].each{ |cc|
      arrCC = cc.split(/-/)
      if( arrCC[0].to_i <= pos2 && pos2 <= arrCC[1].to_i )
        pass_screen = false
      end
    }
    if( unique_read_id.key?(arrSplit[0]) )
      if( unique_read_id[arrSplit[0]] != 2 )
        pass_screen = false
      end
    else
      pass_screen = false
    end
    ##OMIT CONSISTENTLY MAPPED INWARD & OUTWARD FACING READ PAIRS
    if( chrm1 == chrm2 )
      if( arrSplit[8].to_i.abs <= PARAM_INWARD_MAX ) 
        if( (pos1 < pos2) && (ori1 == "+" && ori2 == "-") )
          pass_screen = false
        elsif( (pos1 > pos2) && (ori1 == "-" && ori2 == "+") )
          pass_screen = false
        end
      elsif( arrSplit[8].to_i.abs <= PARAM_OUTWARD_MAX && arrSplit[8].to_i.abs >= PARAM_OUTWARD_MIN )
        if( (pos1 < pos2) && (ori1 =="-" && ori2 == "+") )
          pass_screen = false
        elsif( (pos1 > pos2) && (ori1 == "+" && ori2 == "-") )
          pass_screen = false
        end
      end
    end
    
    if( pass_screen == true )
      pass_counter += 1
      if( pass_counter % 10000 == 0 )
        print "PASS:#{pass_counter/1000}K ( #{((counter/total.to_f)*100).to_i}% ) ... "
      end  
      if(chrm1 < chrm2)
        chr_key = chrm1 + '-' + chrm2
      else
        chr_key = chrm2 + '-' + chrm1
      end
      if(ori1 == ori2)
        ori_key = 'same'
      else
        ori_key = 'diff'
      end
      
      insert_flag = false
      if( clusterHash.key?( chr_key ) )
        if( clusterHash[chr_key].key?( ori_key ) )
          if( clusterHash[chr_key][ori_key].empty? )
            cluster = Cluster.new(line)
            clusterHash[chr_key][ori_key].push(cluster)
          else
            clusterHash[chr_key][ori_key].each_index{ |ii|
              if( clusterHash[chr_key][ori_key][ii].intersect?(line) == true )
                insert_flag = clusterHash[chr_key][ori_key][ii].addBP(line)
                if( insert_flag == true )
                  break
                end
              end
            }
            if( insert_flag == false)
              cluster = Cluster.new(line)
              clusterHash[chr_key][ori_key].push(cluster)
            end
          end
        else
          cluster = Cluster.new(line)
          clusterHash[chr_key][ori_key].push(cluster)
        end
      else
        cluster = Cluster.new(line)
        clusterHash[chr_key][ori_key].push(cluster)
      end
    end
  end
}
print "done\n"
counter = 0
pass_counter = 0
total = 0
cmd.clear
unique_read_id.clear

puts "FINISHED Reading in #{optHash["--BAMFile"]} BAM file and building read clusters."


### CLUSTER QUALITY CONTROL FILTERING
if( optHash.key?("--SV_QC_Filtering") )
  sv_filtering = true
else
  sv_filtering = false
end

if( sv_filtering == true )
  qc_filter_count = 0
  clusterHash.each_key{ |kk|
    clusterHash[kk].each_key{ |oo|
      clusterHash[kk][oo].each_index{ |ii|
        if( clusterHash[kk][oo][ii] != nil && clusterHash[kk][oo][ii].count >= MIN_CLUSTER_COUNT )
          clusterHash[kk][oo].each_index{ |jj|
            if( clusterHash[kk][oo][jj] != nil && jj != ii && clusterHash[kk][oo][jj].count >= MIN_CLUSTER_COUNT )
              if( clusterHash[kk][oo][ii].cluster_intersect?( clusterHash[kk][oo][jj] ) == true )
                if( clusterHash[kk][oo][ii].count > clusterHash[kk][oo][jj].count )
                  clusterHash[kk][oo][jj].qc = 0
                  qc_filter_count += 1
                elsif( clusterHash[kk][oo][ii].count < clusterHash[kk][oo][jj].count )
                  clusterHash[kk][oo][ii].qc = 0
                  qc_filter_count += 1
                elsif( clusterHash[kk][oo][ii].count == clusterHash[kk][oo][jj].count )
                  if( clusterHash[kk][oo][ii].size > clusterHash[kk][oo][jj].size )
                    clusterHash[kk][oo][jj].qc = 0
                    qc_filter_count += 1
                  elsif( clusterHash[kk][oo][ii].size < clusterHash[kk][oo][jj].size )
                    clusterHash[kk][oo][ii].qc = 0
                    qc_filter_count += 1
                  elsif( clusterHash[kk][oo][ii].size == clusterHash[kk][oo][jj].size )
                    clusterHash[kk][oo][ii].qc = 0
                    clusterHash[kk][oo][jj].qc = 0
                    qc_filter_count += 2
                  end
                end
                break
              end
            end
          }
        end
      }
    }
  }
  puts "QUALITY CONTROL FILTERING FLAGGED #{qc_filter_count} CLUSTERS"
end

### MERGING INWARD and OUTWARD FACING READ CLUSTERS
puts "Merging Inward and Outward Facing READ CLUSTERS..."

clusterHash.each_key{ |kk|
  clusterHash[kk].each_key{ |oo|
    puts "Processing #{kk}-#{oo}_ori clusters..."
    clusterHash[kk][oo].each_index{ |ii|
      if( clusterHash[kk][oo][ii] != nil && clusterHash[kk][oo][ii].merge == 0 && clusterHash[kk][oo][ii].count >= MIN_CLUSTER_COUNT && clusterHash[kk][oo][ii].qc == 1 && clusterHash[kk][oo][ii].size > (PARAM_INWARD_MAX * SPAN).to_i )
        outward_read_cluster_size1 = clusterHash[kk][oo][ii].pos1_max - clusterHash[kk][oo][ii].pos1_min + 1
        outward_read_cluster_size2 = clusterHash[kk][oo][ii].pos2_max - clusterHash[kk][oo][ii].pos2_min + 1
        if( outward_read_cluster_size1 >= (seq_length * MIN_CLUSTER_COUNT) && outward_read_cluster_size2 >= (seq_length * MIN_CLUSTER_COUNT) )
          clusterHash[kk][oo].each_index{ |jj|
            if ( clusterHash[kk][oo][jj] != nil && jj != ii && clusterHash[kk][oo][jj].merge == 0 && clusterHash[kk][oo][jj].count >= MIN_CLUSTER_COUNT && clusterHash[kk][oo][jj].qc == 1 && clusterHash[kk][oo][jj].size < (PARAM_INWARD_MAX * SPAN).to_i )
              inward_read_cluster_size1 = clusterHash[kk][oo][jj].pos1_max - clusterHash[kk][oo][jj].pos1_min + 1
              inward_read_cluster_size2 = clusterHash[kk][oo][jj].pos2_max - clusterHash[kk][oo][jj].pos2_min + 1
              if( inward_read_cluster_size1 >= (seq_length * MIN_CLUSTER_COUNT) && inward_read_cluster_size2 >= (seq_length * MIN_CLUSTER_COUNT) )
                if( clusterHash[kk][oo][ii].inward_outward_intersect?(clusterHash[kk][oo][jj]) == true )
                  if( clusterHash[kk][oo][ii].inward_outward_merge(clusterHash[kk][oo][jj]) )
                    puts "Merged INWARD/OUTWARD Clusters = [chrm#{clusterHash[kk][oo][ii].chrm1}(#{clusterHash[kk][oo][ii].ori1}):#{clusterHash[kk][oo][ii].pos1_min}-#{clusterHash[kk][oo][ii].pos1_max} ; chrm#{clusterHash[kk][oo][ii].chrm2}(#{clusterHash[kk][oo][ii].ori2}):#{clusterHash[kk][oo][ii].pos2_min}-#{clusterHash[kk][oo][ii].pos2_max}] , [chrm#{clusterHash[kk][oo][jj].chrm1}(#{clusterHash[kk][oo][jj].ori1}):#{clusterHash[kk][oo][jj].pos1_min}-#{clusterHash[kk][oo][jj].pos1_max} ; chrm#{clusterHash[kk][oo][jj].chrm2}(#{clusterHash[kk][oo][jj].ori2}):#{clusterHash[kk][oo][jj].pos2_min}-#{clusterHash[kk][oo][jj].pos2_max}]"
                    clusterHash[kk][oo][jj] = nil
                    break
                  end
                end
              end
            end
          }
        end
      end
    }
  }
}


### MERGING INVERSION ANNOTATIONS
puts "Combining same Inversion annotations..."

clusterHash.each_key{ |kk|
  chrm_array = kk.split("-")
  if( chrm_array.first == chrm_array.last )
    puts "Processing #{kk}-same_ori INV clusters..."
    clusterHash[kk]["same"].each_index{ |ii|
      if( clusterHash[kk]["same"][ii] != nil && clusterHash[kk]["same"][ii].count >= MIN_CLUSTER_COUNT && clusterHash[kk]["same"][ii].qc == 1 && clusterHash[kk]["same"][ii].size > (PARAM_INWARD_MAX * SPAN).to_i )
        ii_read_cluster_size1 = clusterHash[kk]["same"][ii].pos1_max - clusterHash[kk]["same"][ii].pos1_min + 1
        ii_read_cluster_size2 = clusterHash[kk]["same"][ii].pos2_max - clusterHash[kk]["same"][ii].pos2_min + 1
        if( ii_read_cluster_size1 >= (seq_length * MIN_CLUSTER_COUNT) && ii_read_cluster_size2 >= (seq_length * MIN_CLUSTER_COUNT) )
          if( clusterHash[kk]["same"][ii].getType == "INV" && clusterHash[kk]["same"][ii].inv_merge == 0 )
            clusterHash[kk]["same"].each_index{ |jj|
              if (clusterHash[kk]["same"][jj] != nil && jj != ii && clusterHash[kk]["same"][jj].count >= MIN_CLUSTER_COUNT && clusterHash[kk]["same"][jj].qc == 1 && clusterHash[kk]["same"][jj].size > (PARAM_INWARD_MAX * SPAN).to_i )
                jj_read_cluster_size1 = clusterHash[kk]["same"][jj].pos1_max - clusterHash[kk]["same"][jj].pos1_min + 1
                jj_read_cluster_size2 = clusterHash[kk]["same"][jj].pos2_max - clusterHash[kk]["same"][jj].pos2_min + 1
                if( jj_read_cluster_size1 >= (seq_length * MIN_CLUSTER_COUNT) && jj_read_cluster_size2 >= (seq_length * MIN_CLUSTER_COUNT) )
                  if( clusterHash[kk]["same"][jj].getType == "INV" && clusterHash[kk]["same"][jj].inv_merge == 0 )
                    if( clusterHash[kk]["same"][ii].inverse_intersect?(clusterHash[kk]["same"][jj],kk,"same",ii,jj) == true )
                      puts "Combined same INVERSION annotations = [chrm#{clusterHash[kk]["same"][ii].chrm1}(#{clusterHash[kk]["same"][ii].ori1}):#{clusterHash[kk]["same"][ii].pos1_min}-#{clusterHash[kk]["same"][ii].pos1_max} ; chrm#{clusterHash[kk]["same"][ii].chrm2}(#{clusterHash[kk]["same"][ii].ori2}):#{clusterHash[kk]["same"][ii].pos2_min}-#{clusterHash[kk]["same"][ii].pos2_max}] , [chrm#{clusterHash[kk]["same"][jj].chrm1}(#{clusterHash[kk]["same"][jj].ori1}):#{clusterHash[kk]["same"][jj].pos1_min}-#{clusterHash[kk]["same"][jj].pos1_max} ; chrm#{clusterHash[kk]["same"][jj].chrm2}(#{clusterHash[kk]["same"][jj].ori2}):#{clusterHash[kk]["same"][jj].pos2_min}-#{clusterHash[kk]["same"][jj].pos2_max}]"
                      break
                    end
                  end
                end
              end
            }
          end
        end
      end
    }
  end
}


puts "Writing out Chromosomal Aberrations..."

svp_out = base_name + ".svp"
svp_out_file = File.open(svp_out, "w")
bed_out = base_name + ".bed"
bed_out_file = File.open(bed_out, "w")
bp_out = base_name + ".lff"
bp_out_file = File.open(bp_out, "w")
circos_link = base_name + ".circos.link.txt"
circos_link_file = File.open(circos_link, "w")

sv_count = 0
sv_name_str = nil
if( optHash.key?("--SVname") )
  sv_name_str = optHash["--SVname"]
else
  sv_name_str = "SV"
end

svp_out_file.puts "##program=SVachra-v1.9"
svp_out_file.puts "##abbrev=#{sv_name_str}"
svp_out_file.puts "##source=#{optHash["--BAMFile"]}"
svp_out_file.puts "##META:"
svp_out_file.puts "##META: all SVachra INFO lines will start with the letters #{sv_name_str}"
svp_out_file.puts "##INFO=<ID=\"SVPNAME=#{sv_name_str}00\",Type=String,Description=\"unique name given to this event\">"
svp_out_file.puts "##INFO=<ID=\"#{sv_name_str}TY\",Type=String,Description=\"type\">"
svp_out_file.puts "##INFO=<ID=\"#{sv_name_str}O1\",Type=String,Description=\"orientation 1\">"
svp_out_file.puts "##INFO=<ID=\"#{sv_name_str}O2\",Type=String,Description=\"orientation 2\">"
svp_out_file.puts "##INFO=<ID=\"#{sv_name_str}NR\",Type=Integer,Description=\"num reads\">"
svp_out_file.puts "##INFO=<ID=\"#{sv_name_str}MG\",Type=Boolean,Description=\"Merged Outward and Inward facing reads to call SV (TRUE=1), else only Outward facing reads (FALSE=0)\">"
svp_out_file.puts "##INFO=<ID=\"#{sv_name_str}CTX\",Type=String,Description=\"Coordinates for CTX's pair (chr:pos format)\">"
svp_out_file.puts "#CHROM\tOUTERSTART\tSTART\tINNERSTART\tINNEREND\tEND\tOUTEREND\tTYPE\tSIZE\tREFSEQ\tVARSEQ\tINFO"

clusterHash.each_key{ |kk|
  clusterHash[kk].each_key{ |oo|
    clusterHash[kk][oo].each_index{ |ii|
      if( clusterHash[kk][oo][ii] != nil && clusterHash[kk][oo][ii].count >= MIN_CLUSTER_COUNT && clusterHash[kk][oo][ii].qc == 1 && clusterHash[kk][oo][ii].size > (PARAM_INWARD_MAX * SPAN).to_i )
        read_cluster_size1 = clusterHash[kk][oo][ii].pos1_max - clusterHash[kk][oo][ii].pos1_min + 1
        read_cluster_size2 = clusterHash[kk][oo][ii].pos2_max - clusterHash[kk][oo][ii].pos2_min + 1
        if( read_cluster_size1 >= (seq_length * MIN_CLUSTER_COUNT) && read_cluster_size2 >= (seq_length * MIN_CLUSTER_COUNT) )

          if( clusterHash[kk][oo][ii].getType == "INS" || clusterHash[kk][oo][ii].getType == "DEL" )
            mid1 = (clusterHash[kk][oo][ii].pos1_min + clusterHash[kk][oo][ii].pos1_max)/2
            mid2 = (clusterHash[kk][oo][ii].pos2_min + clusterHash[kk][oo][ii].pos2_max)/2
            sv_size = 0
            avg_frag = 0
            sum_frag = 0

            clusterHash[kk][oo][ii].indel.each{ |frag|
              sum_frag += frag
            }
            avg_frag = (sum_frag / clusterHash[kk][oo][ii].indel.length).to_i
            sv_size = avg_frag.abs
            if( mid1 < mid2 )
              anno_size = clusterHash[kk][oo][ii].pos2_min - clusterHash[kk][oo][ii].pos1_max
              if(anno_size > 0  && sv_size > WINDOW)
                sv_count += 1
                svp_out_file.puts "#{clusterHash[kk][oo][ii].chrm1}\t#{clusterHash[kk][oo][ii].pos1_max}\t.\t.\t.\t.\t#{clusterHash[kk][oo][ii].pos2_min}\t#{clusterHash[kk][oo][ii].getType}\t#{sv_size}\t.\t.\tSVPNAME=#{sv_name_str}#{sv_count};#{sv_name_str}TY=#{clusterHash[kk][oo][ii].getType};#{sv_name_str}O1=#{clusterHash[kk][oo][ii].getOri(1)};#{sv_name_str}O2=#{clusterHash[kk][oo][ii].getOri(2)};#{sv_name_str}NR=#{clusterHash[kk][oo][ii].count};#{sv_name_str}MG=#{clusterHash[kk][oo][ii].merge}"
                bed_out_file.puts "#{clusterHash[kk][oo][ii].chrm1}\t#{clusterHash[kk][oo][ii].pos1_max}\t#{clusterHash[kk][oo][ii].pos2_min}\t#{sv_name_str}#{sv_count}\t#{clusterHash[kk][oo][ii].count}"
                bp_out_file.puts "SVachra\t#{sv_name_str}#{sv_count}\tSV\t#{clusterHash[kk][oo][ii].getType}\tchr#{clusterHash[kk][oo][ii].chrm1}\t#{clusterHash[kk][oo][ii].pos1_min}\t#{clusterHash[kk][oo][ii].pos1_max}\t#{clusterHash[kk][oo][ii].getOri(1)}\t.\t#{clusterHash[kk][oo][ii].count}"
                bp_out_file.puts "SVachra\t#{sv_name_str}#{sv_count}\tSV\t#{clusterHash[kk][oo][ii].getType}\tchr#{clusterHash[kk][oo][ii].chrm2}\t#{clusterHash[kk][oo][ii].pos2_min}\t#{clusterHash[kk][oo][ii].pos2_max}\t#{clusterHash[kk][oo][ii].getOri(2)}\t.\t#{clusterHash[kk][oo][ii].count}"
                circos_link_file.puts "#{sv_name_str}#{sv_count}\ths#{clusterHash[kk][oo][ii].chrm1}\t#{clusterHash[kk][oo][ii].pos1_min}\t#{clusterHash[kk][oo][ii].pos1_max}\tcolor=dgreen"
                circos_link_file.puts "#{sv_name_str}#{sv_count}\ths#{clusterHash[kk][oo][ii].chrm2}\t#{clusterHash[kk][oo][ii].pos2_min}\t#{clusterHash[kk][oo][ii].pos2_max}\tcolor=dgreen"
              end
            elsif( mid2 < mid1 )
              anno_size = clusterHash[kk][oo][ii].pos1_min - clusterHash[kk][oo][ii].pos2_max
              if(anno_size > 0 && sv_size > WINDOW)
                sv_count += 1
                svp_out_file.puts "#{clusterHash[kk][oo][ii].chrm2}\t#{clusterHash[kk][oo][ii].pos2_max}\t.\t.\t.\t.\t#{clusterHash[kk][oo][ii].pos1_min}\t#{clusterHash[kk][oo][ii].getType}\t#{sv_size}\t.\t.\tSVPNAME=#{sv_name_str}#{sv_count};#{sv_name_str}TY=#{clusterHash[kk][oo][ii].getType};#{sv_name_str}O1=#{clusterHash[kk][oo][ii].getOri(2)};#{sv_name_str}O2=#{clusterHash[kk][oo][ii].getOri(1)};#{sv_name_str}NR=#{clusterHash[kk][oo][ii].count};#{sv_name_str}MG=#{clusterHash[kk][oo][ii].merge}"
                bed_out_file.puts "#{clusterHash[kk][oo][ii].chrm2}\t#{clusterHash[kk][oo][ii].pos2_max}\t#{clusterHash[kk][oo][ii].pos1_min}\t#{sv_name_str}#{sv_count}\t#{clusterHash[kk][oo][ii].count}"
                bp_out_file.puts "SVachra\t#{sv_name_str}#{sv_count}\tSV\t#{clusterHash[kk][oo][ii].getType}\tchr#{clusterHash[kk][oo][ii].chrm1}\t#{clusterHash[kk][oo][ii].pos1_min}\t#{clusterHash[kk][oo][ii].pos1_max}\t#{clusterHash[kk][oo][ii].getOri(1)}\t.\t#{clusterHash[kk][oo][ii].count}"
                bp_out_file.puts "SVachra\t#{sv_name_str}#{sv_count}\tSV\t#{clusterHash[kk][oo][ii].getType}\tchr#{clusterHash[kk][oo][ii].chrm2}\t#{clusterHash[kk][oo][ii].pos2_min}\t#{clusterHash[kk][oo][ii].pos2_max}\t#{clusterHash[kk][oo][ii].getOri(2)}\t.\t#{clusterHash[kk][oo][ii].count}"
                circos_link_file.puts "#{sv_name_str}#{sv_count}\ths#{clusterHash[kk][oo][ii].chrm2}\t#{clusterHash[kk][oo][ii].pos2_min}\t#{clusterHash[kk][oo][ii].pos2_max}\tcolor=dgreen"
                circos_link_file.puts "#{sv_name_str}#{sv_count}\ths#{clusterHash[kk][oo][ii].chrm1}\t#{clusterHash[kk][oo][ii].pos1_min}\t#{clusterHash[kk][oo][ii].pos1_max}\tcolor=dgreen"
              end
            end
            
          elsif( clusterHash[kk][oo][ii].getType == "ITX" )
            mid1 = (clusterHash[kk][oo][ii].pos1_min + clusterHash[kk][oo][ii].pos1_max)/2
            mid2 = (clusterHash[kk][oo][ii].pos2_min + clusterHash[kk][oo][ii].pos2_max)/2
            if( mid1 < mid2 )
              anno_size = clusterHash[kk][oo][ii].pos2_max - clusterHash[kk][oo][ii].pos1_max
              if(anno_size > 0)
                sv_count += 1
                svp_out_file.puts "#{clusterHash[kk][oo][ii].chrm1}\t#{clusterHash[kk][oo][ii].pos1_max}\t.\t.\t#{clusterHash[kk][oo][ii].pos2_max}\t.\t.\tMIS\t#{anno_size}\t.\t.\tSVPNAME=#{sv_name_str}#{sv_count};#{sv_name_str}TY=#{clusterHash[kk][oo][ii].getType};#{sv_name_str}O1=#{clusterHash[kk][oo][ii].getOri(1)};#{sv_name_str}O2=#{clusterHash[kk][oo][ii].getOri(2)};#{sv_name_str}NR=#{clusterHash[kk][oo][ii].count};#{sv_name_str}MG=#{clusterHash[kk][oo][ii].merge}"
                bed_out_file.puts "#{clusterHash[kk][oo][ii].chrm1}\t#{clusterHash[kk][oo][ii].pos1_max}\t#{clusterHash[kk][oo][ii].pos2_max}\t#{sv_name_str}#{sv_count}\t#{clusterHash[kk][oo][ii].count}"
                bp_out_file.puts "SVachra\t#{sv_name_str}#{sv_count}\tSV\t#{clusterHash[kk][oo][ii].getType}\tchr#{clusterHash[kk][oo][ii].chrm1}\t#{clusterHash[kk][oo][ii].pos1_min}\t#{clusterHash[kk][oo][ii].pos1_max}\t#{clusterHash[kk][oo][ii].getOri(1)}\t.\t#{clusterHash[kk][oo][ii].count}"
                bp_out_file.puts "SVachra\t#{sv_name_str}#{sv_count}\tSV\t#{clusterHash[kk][oo][ii].getType}\tchr#{clusterHash[kk][oo][ii].chrm2}\t#{clusterHash[kk][oo][ii].pos2_min}\t#{clusterHash[kk][oo][ii].pos2_max}\t#{clusterHash[kk][oo][ii].getOri(2)}\t.\t#{clusterHash[kk][oo][ii].count}"
                circos_link_file.puts "#{sv_name_str}#{sv_count}\ths#{clusterHash[kk][oo][ii].chrm1}\t#{clusterHash[kk][oo][ii].pos1_min}\t#{clusterHash[kk][oo][ii].pos1_max}\tcolor=dgreen"
                circos_link_file.puts "#{sv_name_str}#{sv_count}\ths#{clusterHash[kk][oo][ii].chrm2}\t#{clusterHash[kk][oo][ii].pos2_min}\t#{clusterHash[kk][oo][ii].pos2_max}\tcolor=dgreen"
              end
            elsif( mid2 < mid1 )
              anno_size = clusterHash[kk][oo][ii].pos1_max - clusterHash[kk][oo][ii].pos2_max
              if(anno_size > 0)
                sv_count += 1
                svp_out_file.puts "#{clusterHash[kk][oo][ii].chrm2}\t#{clusterHash[kk][oo][ii].pos2_max}\t.\t.\t#{clusterHash[kk][oo][ii].pos1_max}\t.\t.\tMIS\t#{anno_size}\t.\t.\tSVPNAME=#{sv_name_str}#{sv_count};#{sv_name_str}TY=#{clusterHash[kk][oo][ii].getType};#{sv_name_str}O1=#{clusterHash[kk][oo][ii].getOri(2)};#{sv_name_str}O2=#{clusterHash[kk][oo][ii].getOri(1)};#{sv_name_str}NR=#{clusterHash[kk][oo][ii].count};#{sv_name_str}MG=#{clusterHash[kk][oo][ii].merge}"
                bed_out_file.puts "#{clusterHash[kk][oo][ii].chrm2}\t#{clusterHash[kk][oo][ii].pos2_max}\t#{clusterHash[kk][oo][ii].pos1_max}\t#{sv_name_str}#{sv_count}\t#{clusterHash[kk][oo][ii].count}"
                bp_out_file.puts "SVachra\t#{sv_name_str}#{sv_count}\tSV\t#{clusterHash[kk][oo][ii].getType}\tchr#{clusterHash[kk][oo][ii].chrm1}\t#{clusterHash[kk][oo][ii].pos1_min}\t#{clusterHash[kk][oo][ii].pos1_max}\t#{clusterHash[kk][oo][ii].getOri(1)}\t.\t#{clusterHash[kk][oo][ii].count}"
                bp_out_file.puts "SVachra\t#{sv_name_str}#{sv_count}\tSV\t#{clusterHash[kk][oo][ii].getType}\tchr#{clusterHash[kk][oo][ii].chrm2}\t#{clusterHash[kk][oo][ii].pos2_min}\t#{clusterHash[kk][oo][ii].pos2_max}\t#{clusterHash[kk][oo][ii].getOri(2)}\t.\t#{clusterHash[kk][oo][ii].count}"
                circos_link_file.puts "#{sv_name_str}#{sv_count}\ths#{clusterHash[kk][oo][ii].chrm2}\t#{clusterHash[kk][oo][ii].pos2_min}\t#{clusterHash[kk][oo][ii].pos2_max}\tcolor=dgreen"
                circos_link_file.puts "#{sv_name_str}#{sv_count}\ths#{clusterHash[kk][oo][ii].chrm1}\t#{clusterHash[kk][oo][ii].pos1_min}\t#{clusterHash[kk][oo][ii].pos1_max}\tcolor=dgreen"
              end
            end
            
          elsif( clusterHash[kk][oo][ii].getType == "CTX" )
            sv_count += 1
            svp_out_file.puts "#{clusterHash[kk][oo][ii].chrm1}\t#{clusterHash[kk][oo][ii].pos1_min}\t.\t.\t.\t.\t#{clusterHash[kk][oo][ii].pos1_max}\tMIS\t0\t.\t.\tSVPNAME=#{sv_name_str}#{sv_count};#{sv_name_str}TY=#{clusterHash[kk][oo][ii].getType};#{sv_name_str}O1=#{clusterHash[kk][oo][ii].getOri(1)};#{sv_name_str}O2=#{clusterHash[kk][oo][ii].getOri(2)};#{sv_name_str}NR=#{clusterHash[kk][oo][ii].count};#{sv_name_str}MG=#{clusterHash[kk][oo][ii].merge};#{sv_name_str}CTX=#{clusterHash[kk][oo][ii].chrm2}:#{clusterHash[kk][oo][ii].pos2_min}-#{clusterHash[kk][oo][ii].pos2_max}"
            bed_out_file.puts "#{clusterHash[kk][oo][ii].chrm1}\t#{clusterHash[kk][oo][ii].pos1_min}\t#{clusterHash[kk][oo][ii].pos1_max}\t#{sv_name_str}#{sv_count}\t#{clusterHash[kk][oo][ii].count}"
            bp_out_file.puts "SVachra\t#{sv_name_str}#{sv_count}\tSV\t#{clusterHash[kk][oo][ii].getType}\tchr#{clusterHash[kk][oo][ii].chrm1}\t#{clusterHash[kk][oo][ii].pos1_min}\t#{clusterHash[kk][oo][ii].pos1_max}\t#{clusterHash[kk][oo][ii].getOri(1)}\t.\t#{clusterHash[kk][oo][ii].count}"
            bp_out_file.puts "SVachra\t#{sv_name_str}#{sv_count}\tSV\t#{clusterHash[kk][oo][ii].getType}\tchr#{clusterHash[kk][oo][ii].chrm2}\t#{clusterHash[kk][oo][ii].pos2_min}\t#{clusterHash[kk][oo][ii].pos2_max}\t#{clusterHash[kk][oo][ii].getOri(2)}\t.\t#{clusterHash[kk][oo][ii].count}"
            circos_link_file.puts "#{sv_name_str}#{sv_count}\ths#{clusterHash[kk][oo][ii].chrm1}\t#{clusterHash[kk][oo][ii].pos1_min}\t#{clusterHash[kk][oo][ii].pos1_max}\tcolor=dpurple"
            circos_link_file.puts "#{sv_name_str}#{sv_count}\ths#{clusterHash[kk][oo][ii].chrm2}\t#{clusterHash[kk][oo][ii].pos2_min}\t#{clusterHash[kk][oo][ii].pos2_max}\tcolor=dpurple"
            svp_out_file.puts "#{clusterHash[kk][oo][ii].chrm2}\t#{clusterHash[kk][oo][ii].pos2_min}\t.\t.\t.\t.\t#{clusterHash[kk][oo][ii].pos2_max}\tMIS\t0\t.\t.\tSVPNAME=#{sv_name_str}#{sv_count};#{sv_name_str}TY=#{clusterHash[kk][oo][ii].getType};#{sv_name_str}O1=#{clusterHash[kk][oo][ii].getOri(2)};#{sv_name_str}O2=#{clusterHash[kk][oo][ii].getOri(1)};#{sv_name_str}NR=#{clusterHash[kk][oo][ii].count};#{sv_name_str}MG=#{clusterHash[kk][oo][ii].merge};#{sv_name_str}CTX=#{clusterHash[kk][oo][ii].chrm1}:#{clusterHash[kk][oo][ii].pos1_min}-#{clusterHash[kk][oo][ii].pos1_max}"
            bed_out_file.puts "#{clusterHash[kk][oo][ii].chrm2}\t#{clusterHash[kk][oo][ii].pos2_min}\t#{clusterHash[kk][oo][ii].pos2_max}\t#{sv_name_str}#{sv_count}\t#{clusterHash[kk][oo][ii].count}"
            
          elsif( clusterHash[kk][oo][ii].getType == "INV" )
            if( clusterHash[kk][oo][ii].inv_merge == 0 )
              mid1 = (clusterHash[kk][oo][ii].pos1_min + clusterHash[kk][oo][ii].pos1_max)/2
              mid2 = (clusterHash[kk][oo][ii].pos2_min + clusterHash[kk][oo][ii].pos2_max)/2
              if( mid1 < mid2 )
                if( clusterHash[kk][oo][ii].ori1 == "+" )
                  anno_size = clusterHash[kk][oo][ii].pos2_min - clusterHash[kk][oo][ii].pos1_min
                  if(anno_size > 0)
                    sv_count += 1
                    svp_out_file.puts "#{clusterHash[kk][oo][ii].chrm1}\t.\t.\t#{clusterHash[kk][oo][ii].pos1_min}\t.\t.\t#{clusterHash[kk][oo][ii].pos2_min}\tMIS\t#{anno_size}\t.\t.\tSVPNAME=#{sv_name_str}#{sv_count};#{sv_name_str}TY=#{clusterHash[kk][oo][ii].getType};#{sv_name_str}O1=#{clusterHash[kk][oo][ii].getOri(1)};#{sv_name_str}O2=#{clusterHash[kk][oo][ii].getOri(2)};#{sv_name_str}NR=#{clusterHash[kk][oo][ii].count};#{sv_name_str}MG=#{clusterHash[kk][oo][ii].merge}"
                    bed_out_file.puts "#{clusterHash[kk][oo][ii].chrm1}\t#{clusterHash[kk][oo][ii].pos1_min}\t#{clusterHash[kk][oo][ii].pos2_min}\t#{sv_name_str}#{sv_count}\t#{clusterHash[kk][oo][ii].count}"
                    bp_out_file.puts "SVachra\t#{sv_name_str}#{sv_count}\tSV\t#{clusterHash[kk][oo][ii].getType}\tchr#{clusterHash[kk][oo][ii].chrm1}\t#{clusterHash[kk][oo][ii].pos1_min}\t#{clusterHash[kk][oo][ii].pos1_max}\t#{clusterHash[kk][oo][ii].getOri(1)}\t.\t#{clusterHash[kk][oo][ii].count}"
                    bp_out_file.puts "SVachra\t#{sv_name_str}#{sv_count}\tSV\t#{clusterHash[kk][oo][ii].getType}\tchr#{clusterHash[kk][oo][ii].chrm2}\t#{clusterHash[kk][oo][ii].pos2_min}\t#{clusterHash[kk][oo][ii].pos2_max}\t#{clusterHash[kk][oo][ii].getOri(2)}\t.\t#{clusterHash[kk][oo][ii].count}"
                    circos_link_file.puts "#{sv_name_str}#{sv_count}\ths#{clusterHash[kk][oo][ii].chrm1}\t#{clusterHash[kk][oo][ii].pos1_min}\t#{clusterHash[kk][oo][ii].pos1_max}\tcolor=dgreen"
                    circos_link_file.puts "#{sv_name_str}#{sv_count}\ths#{clusterHash[kk][oo][ii].chrm2}\t#{clusterHash[kk][oo][ii].pos2_min}\t#{clusterHash[kk][oo][ii].pos2_max}\tcolor=dgreen"
                  end
                elsif( clusterHash[kk][oo][ii].ori1 == "-" )
                  anno_size = clusterHash[kk][oo][ii].pos2_max - clusterHash[kk][oo][ii].pos1_max
                  if(anno_size > 0)
                    sv_count += 1
                    svp_out_file.puts "#{clusterHash[kk][oo][ii].chrm1}\t#{clusterHash[kk][oo][ii].pos1_max}\t.\t.\t#{clusterHash[kk][oo][ii].pos2_max}\t.\t.\tMIS\t#{anno_size}\t.\t.\tSVPNAME=#{sv_name_str}#{sv_count};#{sv_name_str}TY=#{clusterHash[kk][oo][ii].getType};#{sv_name_str}O1=#{clusterHash[kk][oo][ii].getOri(1)};#{sv_name_str}O2=#{clusterHash[kk][oo][ii].getOri(2)};#{sv_name_str}NR=#{clusterHash[kk][oo][ii].count};#{sv_name_str}MG=#{clusterHash[kk][oo][ii].merge}"
                    bed_out_file.puts "#{clusterHash[kk][oo][ii].chrm1}\t#{clusterHash[kk][oo][ii].pos1_max}\t#{clusterHash[kk][oo][ii].pos2_max}\t#{sv_name_str}#{sv_count}\t#{clusterHash[kk][oo][ii].count}"
                    bp_out_file.puts "SVachra\t#{sv_name_str}#{sv_count}\tSV\t#{clusterHash[kk][oo][ii].getType}\tchr#{clusterHash[kk][oo][ii].chrm1}\t#{clusterHash[kk][oo][ii].pos1_min}\t#{clusterHash[kk][oo][ii].pos1_max}\t#{clusterHash[kk][oo][ii].getOri(1)}\t.\t#{clusterHash[kk][oo][ii].count}"
                    bp_out_file.puts "SVachra\t#{sv_name_str}#{sv_count}\tSV\t#{clusterHash[kk][oo][ii].getType}\tchr#{clusterHash[kk][oo][ii].chrm2}\t#{clusterHash[kk][oo][ii].pos2_min}\t#{clusterHash[kk][oo][ii].pos2_max}\t#{clusterHash[kk][oo][ii].getOri(2)}\t.\t#{clusterHash[kk][oo][ii].count}"
                    circos_link_file.puts "#{sv_name_str}#{sv_count}\ths#{clusterHash[kk][oo][ii].chrm1}\t#{clusterHash[kk][oo][ii].pos1_min}\t#{clusterHash[kk][oo][ii].pos1_max}\tcolor=dgreen"
                    circos_link_file.puts "#{sv_name_str}#{sv_count}\ths#{clusterHash[kk][oo][ii].chrm2}\t#{clusterHash[kk][oo][ii].pos2_min}\t#{clusterHash[kk][oo][ii].pos2_max}\tcolor=dgreen"
                  end
                end              
              elsif( mid2 < mid1 )
                if( clusterHash[kk][oo][ii].ori2 == "+" )
                  anno_size = clusterHash[kk][oo][ii].pos1_min - clusterHash[kk][oo][ii].pos2_min
                  if(anno_size > 0)
                    sv_count += 1
                    svp_out_file.puts "#{clusterHash[kk][oo][ii].chrm2}\t.\t.\t#{clusterHash[kk][oo][ii].pos2_min}\t.\t.\t#{clusterHash[kk][oo][ii].pos1_min}\tMIS\t#{anno_size}\t.\t.\tSVPNAME=#{sv_name_str}#{sv_count};#{sv_name_str}TY=#{clusterHash[kk][oo][ii].getType};#{sv_name_str}O1=#{clusterHash[kk][oo][ii].getOri(2)};#{sv_name_str}O2=#{clusterHash[kk][oo][ii].getOri(1)};#{sv_name_str}NR=#{clusterHash[kk][oo][ii].count};#{sv_name_str}MG=#{clusterHash[kk][oo][ii].merge}"
                    bed_out_file.puts "#{clusterHash[kk][oo][ii].chrm2}\t#{clusterHash[kk][oo][ii].pos2_min}\t#{clusterHash[kk][oo][ii].pos1_min}\t#{sv_name_str}#{sv_count}\t#{clusterHash[kk][oo][ii].count}"
                    bp_out_file.puts "SVachra\t#{sv_name_str}#{sv_count}\tSV\t#{clusterHash[kk][oo][ii].getType}\tchr#{clusterHash[kk][oo][ii].chrm1}\t#{clusterHash[kk][oo][ii].pos1_min}\t#{clusterHash[kk][oo][ii].pos1_max}\t#{clusterHash[kk][oo][ii].getOri(1)}\t.\t#{clusterHash[kk][oo][ii].count}"
                    bp_out_file.puts "SVachra\t#{sv_name_str}#{sv_count}\tSV\t#{clusterHash[kk][oo][ii].getType}\tchr#{clusterHash[kk][oo][ii].chrm2}\t#{clusterHash[kk][oo][ii].pos2_min}\t#{clusterHash[kk][oo][ii].pos2_max}\t#{clusterHash[kk][oo][ii].getOri(2)}\t.\t#{clusterHash[kk][oo][ii].count}"
                    circos_link_file.puts "#{sv_name_str}#{sv_count}\ths#{clusterHash[kk][oo][ii].chrm2}\t#{clusterHash[kk][oo][ii].pos2_min}\t#{clusterHash[kk][oo][ii].pos2_max}\tcolor=dgreen"
                    circos_link_file.puts "#{sv_name_str}#{sv_count}\ths#{clusterHash[kk][oo][ii].chrm1}\t#{clusterHash[kk][oo][ii].pos1_min}\t#{clusterHash[kk][oo][ii].pos1_max}\tcolor=dgreen"
                  end
                elsif( clusterHash[kk][oo][ii].ori2 == "-" )
                  anno_size = clusterHash[kk][oo][ii].pos1_max - clusterHash[kk][oo][ii].pos2_max
                  if(anno_size > 0)
                    sv_count += 1
                    svp_out_file.puts "#{clusterHash[kk][oo][ii].chrm2}\t#{clusterHash[kk][oo][ii].pos2_max}\t.\t.\t#{clusterHash[kk][oo][ii].pos1_max}\t.\t.\tMIS\t#{anno_size}\t.\t.\tSVPNAME=#{sv_name_str}#{sv_count};#{sv_name_str}TY=#{clusterHash[kk][oo][ii].getType};#{sv_name_str}O1=#{clusterHash[kk][oo][ii].getOri(2)};#{sv_name_str}O2=#{clusterHash[kk][oo][ii].getOri(1)};#{sv_name_str}NR=#{clusterHash[kk][oo][ii].count};#{sv_name_str}MG=#{clusterHash[kk][oo][ii].merge}"
                    bed_out_file.puts "#{clusterHash[kk][oo][ii].chrm2}\t#{clusterHash[kk][oo][ii].pos2_max}\t#{clusterHash[kk][oo][ii].pos1_max}\t#{sv_name_str}#{sv_count}\t#{clusterHash[kk][oo][ii].count}"
                    bp_out_file.puts "SVachra\t#{sv_name_str}#{sv_count}\tSV\t#{clusterHash[kk][oo][ii].getType}\tchr#{clusterHash[kk][oo][ii].chrm1}\t#{clusterHash[kk][oo][ii].pos1_min}\t#{clusterHash[kk][oo][ii].pos1_max}\t#{clusterHash[kk][oo][ii].getOri(1)}\t.\t#{clusterHash[kk][oo][ii].count}"
                    bp_out_file.puts "SVachra\t#{sv_name_str}#{sv_count}\tSV\t#{clusterHash[kk][oo][ii].getType}\tchr#{clusterHash[kk][oo][ii].chrm2}\t#{clusterHash[kk][oo][ii].pos2_min}\t#{clusterHash[kk][oo][ii].pos2_max}\t#{clusterHash[kk][oo][ii].getOri(2)}\t.\t#{clusterHash[kk][oo][ii].count}"
                    circos_link_file.puts "#{sv_name_str}#{sv_count}\ths#{clusterHash[kk][oo][ii].chrm2}\t#{clusterHash[kk][oo][ii].pos2_min}\t#{clusterHash[kk][oo][ii].pos2_max}\tcolor=dgreen"
                    circos_link_file.puts "#{sv_name_str}#{sv_count}\ths#{clusterHash[kk][oo][ii].chrm1}\t#{clusterHash[kk][oo][ii].pos1_min}\t#{clusterHash[kk][oo][ii].pos1_max}\tcolor=dgreen"
                  end
                end
              end
            elsif( clusterHash[kk][oo][ii].inv_merge == 1 )
              if( clusterHash[kk][oo][ii].ori1 == "+" )
                current_pos1 = clusterHash[kk][oo][ii].pos1_min
                current_pos2 = clusterHash[kk][oo][ii].pos2_min
              elsif( clusterHash[kk][oo][ii].ori1 == "-" )
                current_pos1 = clusterHash[kk][oo][ii].pos1_max
                current_pos2 = clusterHash[kk][oo][ii].pos2_max
              end
              if( clusterHash[clusterHash[kk][oo][ii].inv_hash["chrm"]][clusterHash[kk][oo][ii].inv_hash["ori"]][clusterHash[kk][oo][ii].inv_hash["index"]] != nil )
                merge_cluster = clusterHash[clusterHash[kk][oo][ii].inv_hash["chrm"]][clusterHash[kk][oo][ii].inv_hash["ori"]][clusterHash[kk][oo][ii].inv_hash["index"]]
                merge_cluster_size1 = merge_cluster.pos1_max - merge_cluster.pos1_min + 1
                merge_cluster_size2 = merge_cluster.pos2_max - merge_cluster.pos2_min + 1
                if( merge_cluster_size1 >= (seq_length * MIN_CLUSTER_COUNT) && merge_cluster_size2 >= (seq_length * MIN_CLUSTER_COUNT) )
                  if( merge_cluster.ori1 == "+" )
                    merge_pos1 = merge_cluster.pos1_min
                    merge_pos2 = merge_cluster.pos2_min
                  elsif( merge_cluster.ori1 == "-" )
                    merge_pos1 = merge_cluster.pos1_max
                    merge_pos2 = merge_cluster.pos2_max
                  end
                  array = [["current", current_pos1, clusterHash[kk][oo][ii].ori1], ["current", current_pos2, clusterHash[kk][oo][ii].ori2], ["merge", merge_pos1, merge_cluster.ori1], ["merge", merge_pos2, merge_cluster.ori2]]
                  sorted = array.sort{ |a,b| a[1] <=> b[1] } 
                  if( sorted[0][0] != sorted[1][0] && sorted[0][2] != sorted[1][2] && sorted[2][0] != sorted[3][0] && sorted[2][2] != sorted[3][2] )
                    anno_size = sorted[3][1] - sorted[0][1]
                    if(anno_size > 0)
                      sv_count += 1
                      if( clusterHash[kk][oo][ii].merge == 1 || merge_cluster.merge == 1 )
                        merge_flag = 1
                      else
                        merge_flag = 0
                      end
                      svp_out_file.puts "#{clusterHash[kk][oo][ii].chrm1}\t#{sorted[0][1]}\t.\t#{sorted[1][1]}\t#{sorted[2][1]}\t.\t#{sorted[3][1]}\tMIS\t#{anno_size}\t.\t.\tSVPNAME=#{sv_name_str}#{sv_count};#{sv_name_str}TY=#{clusterHash[kk][oo][ii].getType};#{sv_name_str}O1=#{clusterHash[kk][oo][ii].getOri(1)};#{sv_name_str}O2=#{clusterHash[kk][oo][ii].getOri(2)};#{sv_name_str}NR=#{clusterHash[kk][oo][ii].count + merge_cluster.count};#{sv_name_str}MG=#{merge_flag}"
                      bed_out_file.puts "#{clusterHash[kk][oo][ii].chrm1}\t#{sorted[0][1]}\t#{sorted[3][1]}\t#{sv_name_str}#{sv_count}\t#{clusterHash[kk][oo][ii].count + merge_cluster.count}"
                      bp_out_file.puts "SVachra\t#{sv_name_str}#{sv_count}.1\tSV\t#{clusterHash[kk][oo][ii].getType}\tchr#{clusterHash[kk][oo][ii].chrm1}\t#{clusterHash[kk][oo][ii].pos1_min}\t#{clusterHash[kk][oo][ii].pos1_max}\t#{clusterHash[kk][oo][ii].getOri(1)}\t.\t#{clusterHash[kk][oo][ii].count}"
                      bp_out_file.puts "SVachra\t#{sv_name_str}#{sv_count}.1\tSV\t#{clusterHash[kk][oo][ii].getType}\tchr#{clusterHash[kk][oo][ii].chrm2}\t#{clusterHash[kk][oo][ii].pos2_min}\t#{clusterHash[kk][oo][ii].pos2_max}\t#{clusterHash[kk][oo][ii].getOri(2)}\t.\t#{clusterHash[kk][oo][ii].count}"
                      bp_out_file.puts "SVachra\t#{sv_name_str}#{sv_count}.2\tSV\t#{clusterHash[clusterHash[kk][oo][ii].inv_hash["chrm"]][clusterHash[kk][oo][ii].inv_hash["ori"]][clusterHash[kk][oo][ii].inv_hash["index"]].getType}\tchr#{clusterHash[clusterHash[kk][oo][ii].inv_hash["chrm"]][clusterHash[kk][oo][ii].inv_hash["ori"]][clusterHash[kk][oo][ii].inv_hash["index"]].chrm1}\t#{clusterHash[clusterHash[kk][oo][ii].inv_hash["chrm"]][clusterHash[kk][oo][ii].inv_hash["ori"]][clusterHash[kk][oo][ii].inv_hash["index"]].pos1_min}\t#{clusterHash[clusterHash[kk][oo][ii].inv_hash["chrm"]][clusterHash[kk][oo][ii].inv_hash["ori"]][clusterHash[kk][oo][ii].inv_hash["index"]].pos1_max}\t#{clusterHash[clusterHash[kk][oo][ii].inv_hash["chrm"]][clusterHash[kk][oo][ii].inv_hash["ori"]][clusterHash[kk][oo][ii].inv_hash["index"]].getOri(1)}\t.\t#{clusterHash[clusterHash[kk][oo][ii].inv_hash["chrm"]][clusterHash[kk][oo][ii].inv_hash["ori"]][clusterHash[kk][oo][ii].inv_hash["index"]].count}"
                      bp_out_file.puts "SVachra\t#{sv_name_str}#{sv_count}.2\tSV\t#{clusterHash[clusterHash[kk][oo][ii].inv_hash["chrm"]][clusterHash[kk][oo][ii].inv_hash["ori"]][clusterHash[kk][oo][ii].inv_hash["index"]].getType}\tchr#{clusterHash[clusterHash[kk][oo][ii].inv_hash["chrm"]][clusterHash[kk][oo][ii].inv_hash["ori"]][clusterHash[kk][oo][ii].inv_hash["index"]].chrm2}\t#{clusterHash[clusterHash[kk][oo][ii].inv_hash["chrm"]][clusterHash[kk][oo][ii].inv_hash["ori"]][clusterHash[kk][oo][ii].inv_hash["index"]].pos2_min}\t#{clusterHash[clusterHash[kk][oo][ii].inv_hash["chrm"]][clusterHash[kk][oo][ii].inv_hash["ori"]][clusterHash[kk][oo][ii].inv_hash["index"]].pos2_max}\t#{clusterHash[clusterHash[kk][oo][ii].inv_hash["chrm"]][clusterHash[kk][oo][ii].inv_hash["ori"]][clusterHash[kk][oo][ii].inv_hash["index"]].getOri(2)}\t.\t#{clusterHash[clusterHash[kk][oo][ii].inv_hash["chrm"]][clusterHash[kk][oo][ii].inv_hash["ori"]][clusterHash[kk][oo][ii].inv_hash["index"]].count}"
                      circos_link_file.puts "#{sv_name_str}#{sv_count}\ths#{clusterHash[kk][oo][ii].chrm1}\t#{sorted[0][1]}\t#{sorted[1][1]}\tcolor=dgreen"
                      circos_link_file.puts "#{sv_name_str}#{sv_count}\ths#{clusterHash[kk][oo][ii].chrm2}\t#{sorted[2][1]}\t#{sorted[3][1]}\tcolor=dgreen"
                      clusterHash[clusterHash[kk][oo][ii].inv_hash["chrm"]][clusterHash[kk][oo][ii].inv_hash["ori"]][clusterHash[kk][oo][ii].inv_hash["index"]] = nil
                    end                  
                  end
                end
              end
            end
          end
        end
      end
    }
  }
}

clusterHash.clear

exit 0


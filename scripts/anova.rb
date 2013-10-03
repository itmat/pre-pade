#!/usr/bin/env ruby
require 'csv'
require 'rsruby'
require 'optparse'
require 'logger'

# R Functions
####################################################

R = RSRuby.instance()
R.eval_R("scaleData <- function(x) {
  k = (x-min(x))/(max(x)-min(x))*100
  k/100 * 90 + 5
}")

R.eval_R("avgData <- function(d) {
  dd <- filter(d,rep(0.33,3))
  dd[1] <- d[1]
  dd[length(dd)] <- d[length(d)]
  dd
}")

R.eval_R("runAnova <- function(info,csv_file) {
  info
}")

####################################################

# Initialize logger
def setup_logger(loglevel)
  logger = Logger.new(STDERR)
  case loglevel
  when "debug"
    logger.level = Logger::DEBUG
  when "warn"
    logger.level = Logger::WARN
  when "info"
    logger.leve = Logger::INFO
  else
    logger.level = Logger::ERROR
  end
  logger
end

def setup_options(args)
  options = {}

  opt_parser = OptionParser.new do |opts|
    opts.banner = "Usage: anova.rb [options]"
    opts.separator ""
    opts.on("-g", "--gtf_file [GTF_FILE]",:REQUIRED,String,  "gtf_file with gene annotation") do |gtf_file|
      options[:gtf_file] = gtf_file
    end

    opts.on("-v", "--[no-]verbose", "Run verbosely") do |v|
      options[:log_level] = "info"
    end

    opts.on("-d", "--debug", "Run in debug mode") do |v|
      options[:log_level] = "debug"
    end
  end

  opt_parser.parse!(args)

  options
end

class Gene

  def initialize(name, chr, starts, ends, length)
    @name = name
    @chr = chr
    @starts = starts
    @ends = ends
    @length = length
  end

  attr_accessor :name, :chr, :starts, :ends, :length

end

class HTSseq

  def initialize(filename)
    @filename = filename
  end

  attr_accessor :filename



end

def calc_length(starts,ends)
  length = 0
  starts.each_with_index do |start,i|
    length += ends[i] - start
  end
  length
end

def get_genes(gtf_file)
  genes = []
  last_name = ""
  starts = ""
  chr = ""
  starts = []
  ends = []
  length = 0

  File.open(gtf_file).each do |line|
    line.chomp!
    fields = line.split("\t")
    next unless fields[2] == "exon"
    name = fields[-1].split("transcript_id \"")[1].split("\"")[0]
    last_name = name if last_name == ""
    chr = fields[0] if chr == ""
    if name != last_name
      length = calc_length(starts,ends)
      gene = Gene.new(last_name,chr,starts,ends,length)
      genes << gene
      starts = []
      ends = []
      chr = fields[0]
      length = 0
      last_name = name
    end
    starts << fields[3].to_i
    ends << fields[4].to_i
    # First collect information

  end
  length = calc_length(starts,ends)
  gene = Gene.new(last_name,chr,starts,ends,length)
  genes << gene
  genes
end

def run(argv)
  options = setup_options(argv)
  logger = setup_logger(options[:log_level])
  logger.debug(options)
  logger.debug(argv)

  genes = get_genes(options[:gtf_file])


  argv.each do |file|

  end
end

#run(ARGV)




#q_values = R.getQvalue(p_values)
#R.scaleData(numbers).map {|e| e.to_i }



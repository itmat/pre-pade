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

R.eval_R("runAnovaTest <- function(info,csv_file) {
  info
}")

R.eval_R("runAnova <- function(numbers,groups) {
  podwt <- data.frame(num=c(numbers),groups=factor(groups))
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
    logger.level = Logger::INFO
  else
    logger.level = Logger::ERROR
  end
  logger
end

def setup_options(args)
  options = {}

  opt_parser = OptionParser.new do |opts|
    opts.banner = "Usage: anova.rb [options] htseq_file*"
    opts.separator ""

    opts.on("-a", "--gtf_file [GTF_FILE]",:REQUIRED,String,  "gtf_file with gene annotation") do |gtf_file|
      options[:gtf_file] = gtf_file
    end

    opts.on("-v", "--[no-]verbose", "Run verbosely") do |v|
      options[:log_level] = "info"
    end

    opts.on("-d", "--debug", "Run in debug mode") do |v|
      options[:log_level] = "debug"
    end

    opts.on("-g", "--groups [GROUPS]",:REQUIRED,
      String,
      "grouping, for example 3,3 for two groups with each 3 replicates") do |g|
      options[:groups] = g
    end


  end

  args = ["-h"] if args.length == 0
  opt_parser.parse!(args)
  raise "Please specify htseq_files" if args.length == 0
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
    # genename => count
    @genes = {}
    @number_of_fragments = 0
    @no_feature = 0
    @ambiguous = 0
    @too_low_aQual = 0
    @not_aligned = 0
    @alignment_not_unique = 0
  end

  attr_accessor :filename, :genes, :number_of_fragments, :no_feature, :ambiguous,
    :too_low_aQual, :not_aligned, :alignment_not_unique

  def read()
    File.open(@filename).each do |line|
      line.chomp!
      fields = line.split("\t")
      case fields[0]
      when "no_feature"
        @no_feature = fields[1].to_i
      when "ambiguous"
        @mbiguous = fields[1].to_i
      when "too_low_aQual"
        @too_low_aQual = fields[1].to_i
      when "not_aligned"
        @not_aligned = fields[1].to_i
      when "alignment_not_unique"
        @alignment_not_unique = fields[1].to_i
      else
        @genes[fields[0]] = fields[1].to_i
        @number_of_fragments += fields[1].to_i
      end
    end
  end


end

def calc_length(starts,ends)
  length = 0
  starts.each_with_index do |start,i|
    length += ends[i] - start
  end
  length
end

def fpkm(num_fragment,length_transcript,number_mio_of_reads)
  ((num_fragment.to_f/(length_transcript.to_f/1000))/number_mio_of_reads.to_f).to_f
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

def all_fpkm(filenames,genes)
  all_fpkm_values = {}
  filenames.each do |file|
    hts_obj = HTSseq.new(file)
    hts_obj.read
    hts_obj.genes.each_pair do |gene_name, count|
      gene_obj = genes.select { |e| e.name == gene_name }[0]
      all_fpkm_values[gene_name] = [] if !all_fpkm_values.has_key?(gene_name)
      all_fpkm_values[gene_name] << fpkm(count,gene_obj.length,hts_obj.number_of_fragments)
    end
  end
  all_fpkm_values
end

def format_groups(groups)
  groups = groups.split(",").map { |e| e.to_i  }
  rep = []
  feature = []
  groups.each_with_index do |e,i|
    rep[i] = e
    feature[i] = "F#{i}"
  end

  [feature, rep]
end

def run(argv)
  options = setup_options(argv)
  logger = setup_logger(options[:log_level])
  logger.debug(options)
  logger.debug(argv)

  genes = get_genes(options[:gtf_file])

  all_fpkm_values = all_fpkm(argv,genes)

  factors, rep = format_groups(options[:groups])


end

#run(ARGV)




#q_values = R.getQvalue(p_values)
#R.scaleData(numbers).map {|e| e.to_i }



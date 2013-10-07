#!/usr/bin/env ruby
require 'csv'
require 'rsruby'
require 'optparse'
require 'logger'


$logger = Logger.new(STDERR)
# R Functions
####################################################

R = RSRuby.instance()
#RSRuby.set_default_mode(RSRuby::CLASS_CONVERSION)


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

R.eval_R("qVals <- function(p) {
  fxp = ecdf(p)
  q = p/fxp(p)
}")

def initialize_anova(groups)
   R.eval_R(<<-RCOMMAND)
    runAnova <- function(numbers) {
    num = numbers
    groups = factor(c(#{groups}))
    podwt <- data.frame(num = numbers,groups = factor(c(#{groups})))
    fitpodwt <- lm(num ~ groups, data=podwt)
    d_anova <- anova(fitpodwt)
    d_anova[1,5]
  }
RCOMMAND
end

#def initialize_anova(groups,val)
#   R.eval_R(<<-RCOMMAND)
#    runAnova <- function(numbers) {
#    podwt <- data.frame(num = numbers,groups = factor(c(#{groups})),val=c(rep(1,#{val})))
#    fitpodwt <- lm(val~num * groups, data=podwt)
#    d_anova <- anova(fitpodwt)
#    d_anova[1,5]
#  }
#RCOMMAND
#end


def run_anova(numbers,groups,val)
  #R = RSRuby.instance()
 
  #  fitpodwt <- lm(num~groups, data=podwt)
  #  d_anova <- anova(fitpodwt)
  #  d_anova[1,5]
  #}")

  p_value = R.runAnova(numbers)
  p_value
end



####################################################

# Initialize logger
def setup_logger(loglevel)
  
  case loglevel
  when "debug"
    $logger.level = Logger::DEBUG
  when "warn"
    $logger.level = Logger::WARN
  when "info"
    $logger.level = Logger::INFO
  else
    $logger.level = Logger::ERROR
  end
  $logger
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
    $logger.info("Working on #{file}")
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
  #run_anova("1,2,1,3,4,4,10,4,12","rep(\"L\",3),rep(\"I\",3),rep(\"K\",3)",9)
  groups = groups.split(",").map { |e| e.to_i  }
  rep = []
  feature = []
  groups.each_with_index do |e,i|
    rep[i] = e
    feature[i] = "rep(\"F#{i}\",#{e})"
  end
  val = groups.length
  feature = feature.join(",")
  #puts "feature: #{feature}"
  rep = rep.join(",")
  [feature, rep, val]
end

def p_values(all_fpkm_values,feature,val)
  initialize_anova(feature)
  all_p_values = {}
  all_fpkm_values.each_pair do |gene_name, nums|
    #puts "gene_name: #{gene_name}"
    #next unless nums.any? { |e| e != 0 }
    #puts "gene_name: #{gene_name} YES"
    $logger.info("Working on #{gene_name}")
    $logger.info("For these numbers #{nums.join}")
    numbers = nums#.join(",")
    #puts "Nums: #{numbers}"
    p_value = run_anova(numbers,feature,val)
    p_value = 0 if p_value == []
    #puts p_value
    $logger.debug("P_value for #{numbers.join(",")} is #{p_value}")
    all_p_values[gene_name] = p_value
  end
  all_p_values
end

def q_values(all_p_values)
  new_p_values = {}
  all_p_values.each_pair do |k, v|
    next if v == 0
    new_p_values[k] = v
  end

  all_q_vals = R.qVals(new_p_values.values)
  all_with_q_values = {}
  i = 0
  new_p_values.each_pair do |key,value|
    all_with_q_values[key] = [value, all_q_vals[i]]
    i += 1
  end
  all_with_q_values
end

def run(argv)
  options = setup_options(argv)
  $logger = setup_logger(options[:log_level])
  $logger.debug(options)
  $logger.debug(argv)

  genes = get_genes(options[:gtf_file])

  all_fpkm_values = all_fpkm(argv,genes)

  feature, rep, val = format_groups(options[:groups])
  all_p_values = p_values(all_fpkm_values,feature,val)

  all_with_q_values = q_values(all_p_values)
  all_with_q_values.each_pair do |key,value|
    puts "#{key}\t#{value.join("\t")}\t#{all_fpkm_values[key].join("\t")}"
  end
end

run(ARGV)




#q_values = R.getQvalue(p_values)
#R.scaleData(numbers).map {|e| e.to_i }



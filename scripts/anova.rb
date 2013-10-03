#!/usr/bin/env ruby
require 'csv'
require 'rsruby'
require 'optparse'
require 'logger'

# Functions
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

R.anova("runANOVA <- function(d) {
  
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
      options[:verbose] = "info"
    end

    opts.on("-d", "--debug", "Run in debug mode") do |v|
      options[:verbose] = "debug"
    end
  end

  opt_parser.parse!(args)
  
  options
end


options = setup_options(ARGV)
logger = setup_logger(options[:verbose])
logger.debug(options)
logger.debug(ARGV)




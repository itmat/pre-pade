#!/usr/bin/env ruby
require 'csv'
require 'rsruby'

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

####################################################


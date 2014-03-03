#!/usr/bin/env perl
#
# Statistics::Discrete
#
# Chiara Orsini, CAIDA, UC San Diego
# chiara@caida.org
#
# Copyright (C) 2014 The Regents of the University of California.
#
# Statistics::Discrete is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Statistics::Discrete is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Statistics::Discrete.  If not, see <http://www.gnu.org/licenses/>.
#

package Statistics::Discrete;


use 5.014002;
use strict;
use warnings;
use List::Util;
use Math::Round;

our $VERSION = '0.01';

use constant NO_BINNING   => 0;
use constant LIN_BINNING  => 1;
use constant LOG_BINNING  => 2;

use constant LOG_BASE => 10;
use constant DEFAULT_BINS_NUMBER => 10;


# Constructor
sub new {  
  my $class = shift;
  my $self = {};
  $self->{"data"} = []; # anonymous array constructor
  # binning default values
  $self->{"bin-type"} = NO_BINNING;
  $self->{"optimal-binning"} = 1;
  # $self->{"bins"} - bins' structure
  # $self->{"num-bins"} - num bins when optimal binning is off
  # on-demand stored stats
  # $self->{"stats"}{"Desc"} - descriptive statistics
  # $self->{"stats"}{"Dist"} - distributions
  # $self->{"stats"}{"Dist"}{"binned"} - binned distributions
  # bind self and class  
  bless $self, $class;
  return $self;
}

sub add_data {
  my $self = shift;
  my @data_to_add = @_;
  push(@{$self->{"data"}}, @data_to_add);
  # data has changed - stats and bins need to be computed again
  delete $self->{"stats"};
  delete $self->{"bins"};
  return;
}

sub add_data_from_file {
  my $self = shift;
  my $filename = shift;
  my @cols;
  my $f; 
  my $line;
  open($f, "<", $filename) or die "Cannot open " .$filename . ": $!";
  while($line = <$f>) {
    chomp($line);
    if($line =~ /^\s*#/) {
      next;
    }
    @cols = split /\s/ , $line;
    if(scalar @cols == 1) {   # assume list of values
      push(@{$self->{"data"}}, $cols[0]);
    }
    if(scalar @cols == 2) {   # assume list of id - value pairs
      push(@{$self->{"data"}}, $cols[1]);
    }
    # TODO manage other cases
  }
  close($f);
  # data has changed - stats and bins need to be computed again
  delete $self->{"stats"};
  delete $self->{"bins"};
  return;
}



################## Descriptive Statistics ##################

# Minimum value
# compute the minimum value, save the statistic, and return it
sub minimum {
  my $self = shift;
  if(!defined($self->{"stats"}{"Desc"}{"min"})) {
    my $count = $self->count();
    if($count > 0) {
      $self->{"stats"}{"Desc"}{"min"} = List::Util::min(@{$self->{"data"}});
    }
    else {
      $self->{"stats"}{"Desc"}{"min"} = 0;
    }
  }
  return $self->{"stats"}{"Desc"}{"min"};
}


# Maximum value
# compute the minimum value, save the statistic, and return it
sub maximum {
  my $self = shift;
  if(!defined($self->{"stats"}{"Desc"}{"max"})) {
    my $count = $self->count();
    if($count > 0) {
      $self->{"stats"}{"Desc"}{"max"} = List::Util::max(@{$self->{"data"}});
    }
    else {
      $self->{"stats"}{"Desc"}{"max"} = 0;
    }
  }
  return $self->{"stats"}{"Desc"}{"max"};
}


# count the number of samples provided
# save the statistic, and return it
sub count {
  my $self = shift;
  if(!defined($self->{"stats"}{"Desc"}{"count"})) {
    $self->{"stats"}{"Desc"}{"count"} = scalar @{$self->{"data"}};
  }
  return $self->{"stats"}{"Desc"}{"count"};
}


# http://en.wikipedia.org/wiki/Expected_value
#  the weighted average of the possible values, 
#  using their probabilities as their weights
sub mean {
  my $self = shift;
  if(!defined($self->{"stats"}{"Desc"}{"mean"})) {
    my $cumul_value = 0;
    my $count = $self->count();
    my $v;
    foreach $v(@{$self->{"data"}}) {
      $cumul_value += $v;
    }
    if($count > 0) {
      $self->{"stats"}{"Desc"}{"mean"} = $cumul_value / $count;
    }
    else {
      $self->{"stats"}{"Desc"}{"mean"} = 0;     
    }
  }
  return $self->{"stats"}{"Desc"}{"mean"};
}


# http://en.wikipedia.org/wiki/Variance
sub variance {
  my $self = shift;
  if(!defined($self->{"stats"}{"Desc"}{"variance"})) {
    my $mean = $self->mean();
    my $count = $self->count();
    my $cumul_value = 0;    
    my $square_mean = 0;
    my $v;
    foreach $v(@{$self->{"data"}}) {
      $cumul_value += $v**2;
    }
    if($count > 0) {
      $square_mean = $cumul_value / $count;
    }
    $self->{"stats"}{"Desc"}{"variance"} = $square_mean - ($mean**2);
  }
  return $self->{"stats"}{"Desc"}{"variance"};
}


# http://en.wikipedia.org/wiki/Standard_deviation
sub standard_deviation {
  my $self = shift;
  if(!defined($self->{"stats"}{"Desc"}{"sdev"})) { 
    my $variance = $self->variance();
    $self->{"stats"}{"Desc"}{"sdev"} = sqrt($variance);
  }
  return $self->{"stats"}{"Desc"}{"sdev"};
}


################## Distributions Statistics ##################

# http://en.wikipedia.org/wiki/Frequency_distribution
#  the frequency distribution is an arrangement of the values
#  that one or more variables take in a sample. Each entry in
#  the table contains the frequency or count of the occurrences
#  of values within a particular group or interval
sub frequency_distribution {
  my $self = shift;
  if($self->{"bin-type"} == NO_BINNING) { 
    # $self->{"stats"}{"Dist"}{"frequency"}
    return $self->_frequency_distribution();
  }
  else {
    # $self->{"stats"}{"Dist"}{"binned"}{"frequency"}
    return $self->_binned_frequency_distribution();
  }
}


# http://en.wikipedia.org/wiki/Probability_mass_function
#  the probability mass function (pmf) is a function that gives
#  the probability that a discrete random variable is exactly 
#  equal to some value.
sub probability_mass_function {
  my $self = shift;
  if($self->{"bin-type"} == NO_BINNING) { 
    # $self->{"stats"}{"Dist"}{"pmf"}
    return $self->_probability_mass_function();
  }
  else {
    # $self->{"stats"}{"Dist"}{"binned"}{"pmf"}
    return $self->_binned_probability_mass_function();
  }
}


# http://en.wikipedia.org/wiki/Empirical_distribution_function
#  the empirical distribution function, or empirical cdf, is the
#  cumulative distribution function associated with the empirical 
#  measure of the sample. This cdf is a step function that jumps up
#  by 1/n at each of the n data points.
sub empirical_distribution_function {
  my $self = shift;
  if($self->{"bin-type"} == NO_BINNING) { 
    # $self->{"stats"}{"Dist"}{"cdf"}
    return $self->_empirical_distribution_function();
  }
  else {
    # $self->{"stats"}{"Dist"}{"binned"}{"cdf"}
    return $self->_binned_empirical_distribution_function();
  }
}


# http://en.wikipedia.org/wiki/CCDF
#  how often the random variable is above a particular level.
sub complementary_cumulative_distribution_function  {
  my $self = shift;
  if($self->{"bin-type"} == NO_BINNING) { 
    # $self->{"stats"}{"Dist"}{"ccdf"}
    return $self->_complementary_cumulative_distribution_function();
  }
  else {
    # $self->{"stats"}{"Dist"}{"binned"}{"ccdf"}
    return $self->_binned_complementary_cumulative_distribution_function();
  }
}



################## Regular Distributions Statistics ##################


sub _frequency_distribution {
  my $self = shift;
  if(!defined($self->{"stats"}{"Dist"}{"frequency"})) {     
    $self->{"stats"}{"Dist"}{"frequency"} = {};
    my $v;
    foreach $v(@{$self->{"data"}}) {
      if(!defined($self->{"stats"}{"Dist"}{"frequency"}->{$v})) {
	$self->{"stats"}{"Dist"}{"frequency"}->{$v} = 1;
      }
      else {
	$self->{"stats"}{"Dist"}{"frequency"}->{$v}++;
      }
    }
  }
  return $self->{"stats"}{"Dist"}{"frequency"};
}


sub _probability_mass_function {
  my $self = shift;
  if(!defined($self->{"stats"}{"Dist"}{"pmf"})) {  
    my $count = $self->count();
    $self->_frequency_distribution();
    my ($val,$freq);
    foreach $val(keys %{$self->{"stats"}{"Dist"}{"frequency"}}) {
      $freq =  $self->{"stats"}{"Dist"}{"frequency"}{$val};
      $self->{"stats"}{"Dist"}{"pmf"}{$val} = $freq / $count;      
    }
  } 
  return $self->{"stats"}{"Dist"}{"pmf"};
}


sub _empirical_distribution_function {
  my $self = shift;
  if(!defined($self->{"stats"}{"Dist"}{"cdf"})) {  
    my $count = $self->count();
    $self->_frequency_distribution();
    my $val;
    my $cumul_freq = 0;
    foreach $val( sort {$a<=>$b}  keys %{$self->{"stats"}{"Dist"}{"frequency"}}) {
      $cumul_freq +=  $self->{"stats"}{"Dist"}{"frequency"}{$val};
      $self->{"stats"}{"Dist"}{"cdf"}{$val} = $cumul_freq / $count;      
    }
  }   
  return $self->{"stats"}{"Dist"}{"cdf"};
}


sub _complementary_cumulative_distribution_function  {
  my $self = shift;
  # warning, this creates self->{"stats"} if it doesn't exist
  if(!defined($self->{"stats"}{"Dist"}{"ccdf"})) {  
    my $count = $self->count();
    $self->frequency_distribution();
    my $val;
    my $cumul_freq = $count;
    foreach $val( sort {$a<=>$b}  keys %{$self->{"stats"}{"Dist"}{"frequency"}}) {
      $cumul_freq -=  $self->{"stats"}{"Dist"}{"frequency"}{$val};
      $self->{"stats"}{"Dist"}{"ccdf"}{$val} = $cumul_freq / $count;      
    }
  }   
  return $self->{"stats"}{"Dist"}{"ccdf"};
}



################## Binned Distributions Statistics ##################


sub _binned_frequency_distribution {
  my $self = shift;
  if(!defined($self->{"stats"}{"Dist"}{"binned"}{"frequency"})) {     
    # 1. get bins
    $self->bins(); # set $self->{"bins"}
    # 2. get regular frequency distribution
    $self->_frequency_distribution(); # set $self->{"stats"}{"Dist"}{"frequency"}
    # 3. compute binned property
    my $binned_frequency = {};
    my $i = 0;  # vals iterator
    my @sorted_vals = (sort {$a<=>$b} keys %{$self->{"stats"}{"Dist"}{"frequency"}});
    my $num_vals = scalar @sorted_vals;
    my $bi = 0; # bins iterator
    my @sorted_bins = (sort {$a<=>$b} keys %{$self->{"bins"}});
    my $num_bins = scalar @sorted_bins;  
    my $freq_sum = 0;
    for(; $bi < $num_bins; $bi++) { # for each bin
      for(; $i < $num_vals; $i++) {
	if($sorted_vals[$i] < $self->{"bins"}->{$sorted_bins[$bi]}{"right"}) {
	  $freq_sum += $self->{"stats"}{"Dist"}{"frequency"}->{$sorted_vals[$i]}
	}
	else {
	  $binned_frequency->{$sorted_bins[$bi]} = $freq_sum;
	  $freq_sum = 0;
	  last;
	}
      }

      # if last element of vals then we save the current freq sum
      # in the current bin
      if($i == $num_vals) {
	$binned_frequency->{$sorted_bins[$bi]} = $freq_sum;
      }
    }
    $self->{"stats"}{"Dist"}{"binned"}{"frequency"} = $binned_frequency;
  }
  return $self->{"stats"}{"Dist"}{"binned"}{"frequency"};
}


sub _binned_probability_mass_function {
  my $self = shift;
  if(!defined($self->{"stats"}{"Dist"}{"binned"}{"pmf"})) {
    my $count = $self->count();
    my $bins = $self->bins(); # set $self->{"bins"}
    my $bfd = $self->_binned_frequency_distribution();
    my ($val,$freq, $interval);
    foreach $val(keys %{$bfd}) {      
      $freq =  $self->{"stats"}{"Dist"}{"binned"}{"frequency"}{$val};
      $interval = $self->{"bins"}{$val}{"right"} - $self->{"bins"}{$val}{"left"};
      if($interval > 0) { 
	$self->{"stats"}{"Dist"}{"binned"}{"pmf"}{$val} = ($freq / $count) / $interval;      
      }
      else {
	$self->{"stats"}{"Dist"}{"binned"}{"pmf"}{$val} = $freq / $count;
      }
    }
  }
  return $self->{"stats"}{"Dist"}{"binned"}{"pmf"};
}


sub _binned_empirical_distribution_function {
  my $self = shift;
  if(!defined($self->{"stats"}{"Dist"}{"binned"}{"cdf"})) {
    my $bpmf = $self->_binned_probability_mass_function();
    my $val;
    my $cumul_prob = 0;
    foreach $val(sort {$a<=>$b} keys %{$bpmf}) {      
      $cumul_prob +=  $bpmf->{$val};
      $self->{"stats"}{"Dist"}{"binned"}{"cdf"}{$val} = $cumul_prob;
    }
  }
  return $self->{"stats"}{"Dist"}{"binned"}{"cdf"};
}


sub _binned_complementary_cumulative_distribution_function {
  my $self = shift;
  if(!defined($self->{"stats"}{"Dist"}{"binned"}{"ccdf"})) {
    my $bpmf = $self->_binned_probability_mass_function();
    my $val;
    my $cumul_prob = 0;
    my $old_cumul_prob = 0;
    foreach $val(sort {$b<=>$a} keys %{$bpmf}) {         # descending sort   
      $old_cumul_prob = $cumul_prob;
      $cumul_prob += $bpmf->{$val};
      $self->{"stats"}{"Dist"}{"binned"}{"ccdf"}{$val} = $old_cumul_prob;
    }
  }
  return $self->{"stats"}{"Dist"}{"binned"}{"ccdf"};
}


################## Binned Related Functions ##################


sub set_binning_type {
  my $self = shift;
  my $bin_type = shift;
  # check bin-type validity
  if($bin_type != NO_BINNING && $bin_type != LIN_BINNING &&
     $bin_type != LOG_BINNING) {
    die "Wrong binning type provided\n";    
  }
  # if bin type is already set
  if($bin_type == $self->{"bin-type"}) {
    return;
  }
  # if bin-type is different
  $self->{"bin-type"} = $bin_type;
  delete $self->{"num-bins"};
  delete $self->{"bins"};
  delete $self->{"stats"}{"Dist"}{"binned"};
}


sub set_optimal_binning {
  my $self = shift;
  if( $self->{"optimal-binning"} == 0) {
    delete $self->{"bins"};
    delete $self->{"num-bins"};
    delete $self->{"stats"}{"Dist"}{"binned"};
  }
  $self->{"optimal-binning"} = 1;
}


sub set_custom_num_bins {
 my $self = shift;
 my $nb = shift;
 if($self->{"bin-type"} == NO_BINNING) {
   die "Set a binning type before setting the num of bins\n";
 }
 if( $nb < 2 ) {
   die "Wrong number of bins provided: " . $nb . "\n";
 } 
 # check if there is a change
 if( $self->{"optimal-binning"} == 0 and
     defined ($self->{"num-bins"}) and
     $self->{"num-bins"} == $nb) {
   # nothing changed -> nothing to do
 }
 else {
   # update binning parameters
   $self->{"optimal-binning"} = 0;
   $self->{"num-bins"} = $nb;
   delete $self->{"stats"}{"Dist"}{"binned"};
   delete $self->{"bins"};
 }
}


sub bins {
  my $self = shift;
  if(!defined($self->{"bins"})) {
    
    # CASE 1: if no_binning is set we do not
    # save the binning structure
    if( $self->{"bin-type"} == NO_BINNING) {
      my $curbins = {};
      my $v;
      foreach $v(@{$self->{"data"}}) {
	my $ref = $v;
	$curbins->{$v}{"left"} = $v;
	$curbins->{$v}{"right"} = $v;
      }
      return $curbins;
    }
    
    # CASE 2: if a custom number of bins is
    # provided, we lin-/log- bin the support 
    # accordingly
    if($self->{"optimal-binning"} == 0) {
      if($self->{"bin-type"} == LIN_BINNING) {
	$self->{"bins"} = $self->compute_lin_bins($self->{"num-bins"});
      }
      if($self->{"bin-type"} == LOG_BINNING) {
	$self->{"bins"} = $self->compute_log_bins($self->{"num-bins"});
      }
    }

    # CASE 3: optimal binning
    else {
      my $res = $self->_optimal_binning();
      if($res != 0) {
	print "Optimal binning was not possible on this sample, ";
	print "using the NO_BINNING mode" . "\n";
	my $curbins = {};
	my $v;
	foreach $v(@{$self->{"data"}}) {
	  my $ref = $v;
	  $curbins->{$v}{"left"} = $v;
	  $curbins->{$v}{"right"} = $v;
	}
	$self->{"bin-type"} = NO_BINNING;
	$self->{"bins"} = $curbins;
      }
    }
  }
  return $self->{"bins"};
}


################## Optimal Binning Related Functions ##################

# Compute the optimal binning
sub _optimal_binning { 
  my $self = shift;
  # LEGEND:
  # b_* => bins
  # f_* => binned frequency distribution
  # d_* => binned density distribution
  # s_* => sum of density distribution values  
  #
  # start from extreme values ( min_num_bins and max_num_bins)
  # and use the bisection method on the interval
  # till the optimal is found (i.e. the binning such that
  # the sum of the density function values is the closest
  # to 1)
  my $min_num_bins = 2;
  my ($b_min, $f_min, $d_min, $s_min) = $self->_bin_attempt($min_num_bins);
  my $max_num_bins = scalar $self->_full_support();
  my ($b_max, $f_max, $d_max, $s_max) = $self->_bin_attempt($max_num_bins);
  # initial check
  if( $s_min > 1 or $s_max < 1) { 
    return -1;
  }
  my $avg_num_bins;
  my ($b_avg, $f_avg, $d_avg, $s_avg);
  while( ($max_num_bins - $min_num_bins) > 1) {
    $avg_num_bins = Math::Round::round(($max_num_bins + $min_num_bins)/2);
    ($b_avg, $f_avg, $d_avg, $s_avg) = $self->_bin_attempt($avg_num_bins);
    if($s_avg < 1) {
      $min_num_bins = $avg_num_bins;
      ($b_min, $f_min, $d_min, $s_min) = ($b_avg, $f_avg, $d_avg, $s_avg);
    }
    else {
      $max_num_bins = $avg_num_bins;
      ($b_max, $f_max, $d_max, $s_max) = ($b_avg, $f_avg, $d_avg, $s_avg);
    }
  }
  # the binning whose density sum is closer to 1 win
  if( (1 - $s_min) < ($s_max-1) ) {
    $self->{"bins"} = $b_min;
    $self->{"Dist"}{"binned"}{"frequency"} = $f_min;
    $self->{"Dist"}{"binned"}{"pmf"} = $d_min;
  }
  else {
    $self->{"bins"} = $b_max;
    $self->{"Dist"}{"binned"}{"frequency"} = $f_max;
    $self->{"Dist"}{"binned"}{"pmf"} = $d_max;
  }
  return 0;
}


# bin the current data into $num_bins intervals
# and return: the bins, the binned frequency distribution
# the binned density distribution and the sum of the 
# density values
sub _bin_attempt {
  my $self = shift;
  my $num_bins = shift;
  # output values
  my ($cur_bins, $cur_freq, $cur_density, $cur_sum);
  $cur_sum = 0;
  # get the bins
  if( $self->{"bin-type"} == LIN_BINNING ) {
    $cur_bins = $self->compute_lin_bins($num_bins);
  }
  if( $self->{"bin-type"} == LOG_BINNING ) {
    $cur_bins  =$self->compute_log_bins($num_bins);
  }
  # assert non binned frequency distribution is computed
  $self->_frequency_distribution(); # set $self->{"stats"}{"Dist"}{"frequency"}
  my $count = $self->count();
  #
  my $s = 0; # support iterator
  my @sorted_support = (sort {$a<=>$b} $self->_full_support());
  my $sup_size = scalar @sorted_support;
  my $bi = 0; # bins iterator
  my @sorted_bins = (sort {$a<=>$b} keys %{$cur_bins});
  my $n_bins = scalar @sorted_bins;  # assert(n_bins == num_bins)
  my $ref = 0;
  my $freq_sum = 0;
  my $interval = 0;
  for(; $bi < $n_bins; $bi++) { # for each bin
    $ref = $sorted_bins[$bi];
    $interval = $cur_bins->{$ref}{"right"} - $cur_bins->{$ref}{"left"};
    for(; $s < $sup_size; $s++) { # for each support value
      # if the current support value is less than the rightmost value in the bin
      if($sorted_support[$s] < $cur_bins->{$ref}{"right"}) {
	$freq_sum += $self->{"stats"}{"Dist"}{"frequency"}->{$sorted_support[$s]};
      }
      # else we have to save the information for the current bin
      else {
	$cur_freq->{$ref} = $freq_sum;
	if($interval > 0) {
	  $cur_density->{$ref} = ($freq_sum / $count) / $interval;
	}
	else {
	  $cur_density->{$ref} =  $freq_sum / $count;
	}
	$cur_sum += $cur_density->{$ref};
	$freq_sum = 0;
	last; # we do not increment s
      }
    }    
    # if last element of the support has been reached 
    # we have to save the current freq sum in the current bin
    if($s == $sup_size) {
      $cur_freq->{$ref} = $freq_sum;
      if($interval > 0) {
	$cur_density->{$ref} = ($freq_sum / $count) / $interval;
      }
      else {
	$cur_density->{$ref} =  $freq_sum / $count;
      }
      $cur_sum += $cur_density->{$ref};	
    }
  }
  #DEBUG print "attempt: " . $num_bins . " sum: " . $cur_sum . "\n";
  return ($cur_bins, $cur_freq, $cur_density, $cur_sum);
}



sub _full_support {
  my $self = shift;
  my $support = {};
  my $v;
  foreach $v(@{$self->{"data"}}) {
    $support->{$v} = 1;
  }
  return (keys %{$support});
}




################## Utilities ##################

# return a structure describing a linear binning
# of the current data ($num_bins provided)
sub compute_lin_bins {
  my $self = shift;
  my $num_bins = shift;
  my $min = $self->minimum();
  my $max = $self->maximum();
  my $linbins = {}; # reference to empty hash
  # last bin includes just the maximum value
  my $delta = abs($max - $min)/($num_bins-1);
  my $i;
  my $ref;
  for($i=0; $i < $num_bins; $i++) {
    $ref = $min + $i * $delta + $delta/2;
    $linbins->{$ref}{"left"} =  $min + $i * $delta;
    $linbins->{$ref}{"right"} = $min + $i * $delta + $delta;
  }
  return $linbins;
}


# return a structure describing a logarithmic binning
# of the current data ($num_bins provided)
# bins are not set
sub compute_log_bins {
  my $self = shift;
  my $num_bins = shift;
  my $min = $self->minimum();
  my $max = $self->maximum();
  my $logbins = {}; # reference to empty hash
  my $min_exp = log($min)/log(LOG_BASE);
  my $max_exp = log($max+1)/log(LOG_BASE);
  my $delta_exp = abs($max_exp - $min_exp)/$num_bins;
  my $cur_exp;
  my $i;
  my $ref;
  for($i=0; $i < $num_bins; $i++) {
    $cur_exp = $min_exp + $i * $delta_exp;
    $ref = LOG_BASE ** ($cur_exp + $delta_exp/2);
    $logbins->{$ref}{"left"} = LOG_BASE ** $cur_exp;
    $logbins->{$ref}{"right"} = LOG_BASE ** ($cur_exp + $delta_exp);    
  }
  return $logbins;
}



1;
__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

Statistics::Discrete - Perl extension for blah blah blah

=head1 SYNOPSIS

  use Statistics::Discrete;
 
  # construct a new class
  my $sd = Statistics::Discrete->new();
  
  # add data using an array
  $sd->add_data((2,5,7,2,1,7,3,3,7,333));
  
  # ---- Descriptive Statistics ---- #

  # Count the number of samples
  my $count = $sd->count(); 
  # Minimum data value
  my $min = $sd->minimum(); 
  # Maximum data value
  my $max = $sd->maximum(); 
  # Mean 
  my $mean = $sd->mean(); 
  # Variance
  my $variance = $sd->variance(); 
  # Standard Deviation
  my $standard_deviation = $sd->standard_deviation(); 


  # ---- Distributions ---- #

  my $fd = $sd->frequency_distribution();
  my $pmf = $sd->probability_mass_function();
  my $cdf = $sd->empirical_distribution_function();
  my $ccdf = $sd->complementary_cumulative_distribution_function();


  # ---- Binning ---- #
  # binning influences the way the distributions
  # are returned
  my $binning_type = NO_BINNING;   # DEFAULT
  $binning_type = LIN_BINNING; 
  $binning_type = LOG_BINNING; 
  $sd->set_binning_type($binning_type);

  $sd->set_optimal_binning();
  my $num_of_bins = 3;
  $sd->set_custom_num_bins($num_of_bins);
  
  # return the bins currently used
  my $cur_bins = $sd->bins();

  # compute linear bins starting from the data
  my $lin_bins = $sd->compute_lin_bins($num_of_bins);
  my $log_bins = $sd->compute_log_bins($num_of_bins);

  # ---- Data Input ---- #

  # add data from file
  $sd->add_data_from_file("./data.txt"); 


=head1 DESCRIPTION

Stub documentation for Statistics::Discrete, created by h2xs. 

TODO: work in progress

=head2 EXPORT

None by default.



=head1 SEE ALSO

Mention other useful documentation such as the documentation of
related modules or operating system documentation (such as man pages
in UNIX), or any relevant external documentation such as RFCs or
standards.

If you have a mailing list set up for your module, mention it here.

If you have a web site set up for your module, mention it here.

=head1 AUTHOR

Chiara Orsini, E<lt>chiara@caida.orgE<gt>

=head1 COPYRIGHT AND LICENSE

Chiara Orsini, CAIDA, UC San Diego
chiara@caida.org

Copyright (C) 2014 The Regents of the University of California.

Statistics::Discrete is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Statistics::Discrete is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Statistics::Discrete.  If not, see <http://www.gnu.org/licenses/>.


=cut

package BioStats;

=head1 NAME

BioStats - set of utility functions for analyzing biological data.

=head1 DESCRIPTION

=cut

=head3 getE

<float> getE()

Returns the numerical constant e.

=cut

sub getE {
  return 2.7182818284;
}

=head3 getPi

<float> getPi()

Returns the numerical constant pi.

=cut

sub getPi {
  return 3.1415926535;
}

=head3 getPoissonDistribution

<\@> getPoissonDistribution(<reference to list>;<int>)

Returns a reference to the Poisson distribution over the first n terms given an array of occurances and an integer specifying n. n will be 10 unless specified.

=item my @list = [1, 4, 5, 12, 20];

=item my $median = BioStats::getPoissonDistribution(\@list,10);

=item print $median,"\n";

=back

=cut

sub getPoissonDistribution {
  my ($data,$length) = @_;
  
  my @result;
  $length = 10 if !defined($length) || $length < 0;
  my $mean = getMean($data);
  my $count = @{$data};
  my $base = $count/(getE()**$mean);
  for (my $i = 1; $i < $length; $i++) {
    $result[$i-1] = $base;
    $base = $base * $mean / $i;
  }
  $result[$length-1] = $base;
  
  return \@result;
}

=head3 getMedian

<float> getMedian(<reference to list>)

Returns the median value from a pre-sorted list of numeric values.

=item my @list = [1, 4, 5, 12, 20];

=item my $median = BioStats::getMedian(\@list);

=item print $median,"\n";

=back

The above example prints the value 5.

=cut

sub getMedian {
  my ($data,$sort) = @_;
  $sort = 0 unless defined($sort);
  
  return getMedialValue(0.5,$data,$sort);
}

=head3 getMedialValue

=item <float> getMedialValue(<float>,<reference to list>,<boolean>)

=back

The function getMedialValue returns a value from a specified location within
a list of values. The position is specified as a decimal number indicating
the proportion of values to ignore before selecting a value.

This is useful for accessing the median value of a list. The value returned
depends on the number of elements in the list. For instance, if you want
the median value of a three element list, the value of the middle element
will be returned. However, if you want the median value of a 4 element
array, the average of the 2nd and 3rd elements will be returned.

Empty lists and lists with only a single element are degenerate cases and
return the value 0 or the value of the only element respectively.

=head2 Examples

    Get the median value from a pre-sorted list
    $median = BioStats::getMedialValue(0.5,\@list);
    
    Get the first quadril from an unsorted list
    $firstQuadril = getBiostats::getMedialValue(0.25,\@list,1);
    
=cut

sub getMedialValue {
  my ($percent,$data,$sort) = @_;
  $sort = 0 unless defined($sort);
  
  # small arrays are special cases
  my $size = @{$data};
  if ($size == 0) {
    return 0;
  } elsif ($size == 1) {
    return ${$data}[0];
  }
  
  # if sorting is specified then sort data numerically in ascending order
  my @sortedData;
  if (defined($sort) && !$sort) {
    @sortedData = @$data;
  } else {
    @sortedData = sort { $a <=> $b;} @{$data};
  }
  
  my $ipos = int($size * $percent);
  my $npos = int(($size - 1) * $percent);
  if ($ipos == $npos) {
    # indicates that the position cooresponds to a single element of
    # the list
    return $sortedData[$ipos];
  } else {
    # position falls between two elements so take average of both
    return (($sortedData[$ipos] + $sortedData[$npos]) / 2);
  }
}

=head3 getJulesCantor

=item <float> getJukesCantor(<float>);

=back

Given that there are only 4 nucleotides to choose from, and assuming
stochasticity, with a random assigment of nucleotides to two taxa you'd
expect them to be 75% different. The Jukes-Cantor adjustment provides a
correction to the observed distance between two sequences under these
assumptions. The equation is :

=item Jukes Cantor = - 3/4 ln(1-(4p/3))

=head2 See Also

http://research.amnh.org/users/siddall/methods/day6.html

=cut

sub getJukesCantor {
  my ($distance) = @_;

  return (-3 * log (1 - (4 * $distance / 3)) / 4);
}

=head3 getMean

=item <float> getMean(<reference to list>);

=back

Get the mean or average value from a list of values.

=cut

sub getMean {
  my ($values) = @_;
  
  my $count = 0;
  my $sum = 0;
  foreach my $value (@{$values}) {
    $count++;
    $sum+= $value;
  }
  my $avg;
  if ($count) {
    $avg = $sum / $count;
  } else {
    
  }
  return ($avg);
}

=head3 getMeanModuloOutliers

=item <float>,<float> getMean(<reference to list>,<float>,<float>,<boolean>);

=back

Get the mean or average value from a list of values.

=cut

sub getMeanStdModuloOutliers {
  my ($values,$min,$max,$sort) = @_;
  
  # if sorting is specified then sort data numerically in ascending order
  my @sortedData;
  if (defined($sort) && !$sort) {
    @sortedData = @$values;
  } else {
    @sortedData = sort { $a <=> $b;} @{$values};
  }
	my $size = @sortedData;
	my @final;
	for (my $i = 0; $i < @sortedData; $i++) {
		my $p = $i/$size;
		push @final,$sortedData[$i] if $p >= $min && $p <= $max;
	}
	
	my $mean = getMean(\@final);
	my $std = getMeanStandardDeviation(\@final,$mean);
  return ($mean,$std);
}

=head3 getMeanStandardDeviation

=item <float> getMeanStandardDeviation(<reference to list>[,<float>]);

=back

Get the standard deviation from the mean given a list of values. Optionally
can provide the mean.

=cut

sub getMeanStandardDeviation {
  my ($values,$mean) = @_;
  
  if (!defined($mean)) {
    $mean = getMean($values);
  }
  
  my $sumOfSquares = 0;
  my $count = 0;
  foreach my $value (@{$values}) {
    $sumOfSquares += ($value - $mean) * ($value - $mean);
    $count++;
  }
  my $std = 0;
  if ($count > 1) {
    $std = sqrt($sumOfSquares/($count - 1));
  }
  
  return ($std);
}

=head3 getVariance

=item <float> getVariance(<reference to list>[,<float>]);

=back

Get the variance from the mean given a list of values. Optionally
can provide the mean.

=cut

sub getVariance {
  my ($values,$mean) = @_;
  
  if (!defined($mean)) {
    $mean = getMean($values);
  }
  
  my $sumOfSquares = 0;
  my $count = 0;
  foreach my $value (@{$values}) {
    $sumOfSquares += ($value - $mean) * ($value - $mean);
    $count++;
  }
	my $var = 0;
  if ($count > 1) {
    $var = $sumOfSquares/$count;
  }
  
  return ($var);
}

=head3 zScoreTransform

=item <reference to list> zScoreTransform(<reference to list>[,<positive int>,<float>,<float>]);

=back

Given a list of values, will compute the z-score.

=cut

sub zScoreTransform {
  my ($values,$rounding,$mean,$std) = @_;
  
	if (defined($rounding)) {
		$rounding = 10**$rounding;
	} else {
		$rounding = 1000;
	}
  if (!defined($mean)) {
    $mean = getMean($values);
  }
  if (!defined($std)) {
    $std = getMeanStandardDeviation($values,$mean);
  }
  
	my @result;
  foreach my $value (@{$values}) {
		my $z = int($rounding*($value - $mean)/$std)/$rounding;
		push @result,$z;
  }
  
  return (\@result);
}

=head3 getMin

=item <float> getMin(<reference to list>);

=back

Return the lowest value from a list of values.

=cut

sub getMin {
  my ($values) = @_;
  
  my $min;
  foreach my $value (@{$values}) {
    $min = $value unless defined($min);
    $min = $value if $value < $min;
  }

  return ($min);
}

=head3 getMax

=item <float> getMax(<reference to list>);

=back

Return the max value from a list of values.

=cut

sub getMax {
  my ($values) = @_;
  
  my $max;
  foreach my $value (@{$values}) {
    $max = $value unless defined($max);
    $max = $value if $value > $max;
  }

  return ($max);
}

1;

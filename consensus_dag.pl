#!/usr/bin/perl

use strict;
use warnings;
use DAG;
use AI::Genetic::Pro;
use Getopt::Long;
use Algorithm::Permute;
use Chart::Gnuplot;

use constant MAX_POPULATION => 200;
use constant INITIAL_SCORE => 10000;

# Input number of dags to be used to get consensus dag
my $num_dags = 3;

# Input number of generations to evolve
my $num_gen = 20;
my $file_prefix = 'dag';
my $all = 0;
my $chart_name = 'chart.gif';
my $random = 260;

GetOptions(
    "ndags=i"  => \$num_dags,
    "gen=i"    => \$num_gen,
    "prefix=s" => \$file_prefix,
    "chart:s"  => \$chart_name,
    "random:i" => \$random,
    "all"      => \$all
) or die "Invalid command line arguments";

# Input all dags
my $dags = [];
for(my $i=0; $i < $num_dags; ++$i) {
    push @$dags, new DAG("$file_prefix$i".'.dat');
}

my @vertices = $dags->[0]->{graph}->vertices;
my $fittest = [];
my $best_score = INITIAL_SCORE;

if(!$all) {
    my $population = fact(scalar @vertices);
    $population = ($population < MAX_POPULATION ? $population : MAX_POPULATION);

    my $ga = AI::Genetic::Pro->new(        
        -fitness         => \&fitness,        # fitness function
        -terminate       => \&terminate,      # terminate function
        -type            => 'combination',    # type of chromosomes
        -population      => $population,      # population
        -crossover       => 0.9,              # probab. of crossover
        -mutation        => 0.01,             # probab. of mutation
        -parents         => 2,                # number  of parents
        -selection       => [ 'Roulette' ],   # selection strategy
        -strategy        => [ 'OX' ],         # crossover strategy
        -cache           => 1,                # cache results
        -history         => 1,                # remember best results
        -preserve        => 3,                # remember the bests
        -variable_length => 0,                # turn variable length OFF
	);

    $ga->init( \@vertices );
    my $seen = {};
    foreach(@{$ga->people}) {
	my $t = join(',', @$_);
	$seen->{$t} = 1;
    }
    my $cnt = scalar(keys %$seen);
    print 'Initial chromosomes = ', $cnt, "\n";
    for(my $i=1; $i <= $num_gen; $i++) {
	#$ga->evolve($num_gen);
	$ga->evolve(1);
	foreach(@{$ga->people}) {
	    my $t = join(',', @$_);
	    if(!exists $seen->{$t}) {
		$seen->{$t} = 1;
		$cnt++;
	    }
	}
	print 'Total population = ', $ga->population, "\n";
	print "After $i generations = ", scalar(keys %$seen), "\n";
    }

    #plot_chart($ga, $chart_name);

    $fittest = $ga->getFittest;
    $best_score = 1/$ga->as_value($fittest);
}
else {
    my $p_iterator = Algorithm::Permute->new ( \@vertices );
    my $total_cost = 0;
    my $i = 1;
    while (my @perm = $p_iterator->next) {
	if($i%100 == 0) {
	    print STDERR "Calculating for order $i- $best_score\n";
	}
	foreach my $dag (@$dags) {
	    my $dag_alpha = method_b2($dag, \@perm);
	    $total_cost += $dag_alpha->get_cost();
	}
	if($total_cost < $best_score) {
	    $best_score = $total_cost;
	    $fittest = \@perm;
	}
	$i++;
    }
}

print STDERR 'Best Score: ', $best_score, "\n";
print @{$fittest}, "\n";

for(my $i=0; $i < $num_dags; ++$i) {
    my $dag_alpha = method_b2($dags->[$i], $fittest);
    print "dag$i = ", $dag_alpha->to_string(), "\n";
}

=head2 plot_chart($ga, $chart_name)
=cut
sub plot_chart {

    my ($ga, $chart_name) = @_;

    my $history = $ga->getHistory();

    my $scores = [];
    for(my $i=0; $i < scalar(@$dags); $i++) {
	my $score = fitness( {}, [ $dags->[$i]{graph}->topological_sort ] );
	my $score_array = [];
	foreach(1..$num_gen) {
	    push @$score_array, $score;
	}
	my $dataset = Chart::Gnuplot::DataSet->new(
	    xdata => [1..$num_gen],
	    ydata => $score_array,
	    title => 'Score for graph '.$i,
	    style => "linespoints",
	    );
	push @$scores, $dataset;
    }
    my $chart = Chart::Gnuplot->new(
        output => $chart_name,
        title  => "Score Improvement",
        xlabel => "generation",
        ylabel => "best score",
    );
    my $data = Chart::Gnuplot::DataSet->new(
        xdata => [1..$num_gen],
        ydata => $history->[0],
        title => 'GA best score',
        style => "lines",
    );
    unshift @$scores, $data;
    push @$scores, get_random_data_score();
    $chart->plot2d(@$scores);
}

=head2 get_random_data_score()
=cut
sub get_random_data_score {

    # using Genetic::Pro to calculate random best score by evolving for only 1 generation
    my $ga = AI::Genetic::Pro->new(        
        -fitness         => \&fitness,        # fitness function
        -terminate       => \&terminate,      # terminate function
        -type            => 'combination',    # type of chromosomes
        -population      => $random,          # population
        -crossover       => 0.9,              # probab. of crossover
        -mutation        => 0.01,             # probab. of mutation
        -parents         => 2,                # number  of parents
        -selection       => [ 'Roulette' ],   # selection strategy
        -strategy        => [ 'OX' ],         # crossover strategy
        -cache           => 0,                # cache results
        -history         => 1,                # remember best results
        -preserve        => 3,                # remember the bests
        -variable_length => 0,                # turn variable length OFF
	);

    $ga->init( \@vertices );
    #$ga->evolve(1);
    my $best_score = $ga->as_value($ga->getFittest);
    my $score_data = [];
    for(1..$num_gen) {
	push @$score_data, $best_score;
    }
    my $data = Chart::Gnuplot::DataSet->new(
        xdata => [1..$num_gen],
        ydata => $score_data,
        title => 'Random best score',
        style => "linespoints",
    );
    return $data;
}

=head2 fact($n)
=cut
sub fact {
  my $n = shift;
  $n == 0 ? 1 : $n*fact($n-1);
}

=head2 fitness($ga, $s_alpha)
=cut
sub fitness {
    my ($ga, $s_alpha) = @_;

    my $total_cost = 0;
    foreach my $dag (@$dags) {
	my $dag_alpha = method_b2($dag, $s_alpha);
	$total_cost += $dag_alpha->get_cost();
    }
    $total_cost = ($total_cost ? 1/$total_cost : $total_cost);

    return $total_cost;
}

=head2 terminate($ga)
=cut
sub terminate {
    my ($ga) = @_;
    my $score = $ga->as_value($ga->getFittest);
    return $score <= 0 ? 1 : 0;
}

=head2 occurs_right($s_alphs, $a, $b)
=cut
sub occurs_right {
    my ($s_alpha, $a, $b) = @_;

    my $ret_val = 0;
    foreach my $node (@$s_alpha) {
        # if 'b' comes before 'a', then 'a' is to the right of 'b'
	if($b == $node) {
	    $ret_val = 1;
	    last;
	}
        # if 'a' comes before 'b', then 'a' is to the left to 'b'
	if($a == $node) {
	    $ret_val = 0;
	    last;
	}
    }
    return $ret_val;
}

sub find
{
    my ($s_beta, $value) = @_;
    my $position = 0;
    for(my $i=0; $i < @$s_beta; $i++) {
	if($s_beta->[$i] == $value) {
	    $position = $i;
	    last;
	}
    }
    return $position;
}

sub swap {
    my ($s_beta, $a, $b) = @_;

    my $a_index = find($s_beta, $a);
    my $b_index = find($s_beta, $b);

    # swap 'a' & 'b' in s_beta
    ( $s_beta->[$a_index], $s_beta->[$b_index] ) = ( $s_beta->[$b_index], $s_beta->[$a_index] );
}

sub construct {
    my($g, $s_alpha) = @_;

    # create empty array ordering s_beta
    my $s_beta = [];

    # create deep copy of current DAG
    my $g_dash = $g->create_copy();

    while( scalar( $g_dash->get_vertices() ) > 0 ) {

	# get rightmost sink node from ordering s_alpha
	my $a = $g_dash->get_rightmost_sink_node($s_alpha);

	# add sink_node as leftmost node in s_beta
	if(defined $a) {
	    unshift @$s_beta, $a;
	}

	my $swapped;
	do {
	    $swapped = 0;

	    # Let b denote right neighbor of a in s_beta
	    my $b = undef;
	    if(@$s_beta > 1) {
		$b = $s_beta->[1];
	    }

	    # if 'b' not empty and 'a' is not parent of 'b' in G and 'a' is right of 'b'
	    # in s_alpha, then interchange 'a' & 'b' in s_beta
	    if( defined $b && !$g->is_parent($b, $a) 
		&& occurs_right($s_alpha, $a, $b) ) {
		swap($s_beta, $a, $b);
		$swapped = 1;
	    }
	} while($swapped);

	$g_dash->remove_node($a);
	#print scalar($g_dash->get_vertices()), "\n";
    }
    return $s_beta;
}

sub equal {
    my ($arr1, $arr2) = @_;

    if( scalar(@$arr1) == scalar(@$arr2) ) {
	my $cnt = scalar(@$arr1);
	for(my $i=0; $i < $cnt; $i++) {
	    if( $arr1->[$i] != $arr2->[$i] ) {
		return 0;
	    }
	}
    }
    else {
	return 0;
    }
    return 1;
}

sub get_rightmost_node {
    my ($s_alpha, $considered) = @_;

    my $right_node = undef;
    foreach my $node (reverse @$s_alpha) {
	if( !$considered->{$node} ) {
	    $right_node = $node;
	    $considered->{$node} = 1;
	    last;
	}
    }
    return $right_node;
}

sub get_right_neighbor {
    my ($s_beta, $y) = @_;
    my $right_neighbor = undef;
    my $position = find($s_beta, $y);
    if( defined $position && $position+1 < scalar(@$s_beta) ) {
	$right_neighbor = $s_beta->[$position+1];
    }
    return $right_neighbor;
}

sub method_b2 {
    my ($dag, $s_alpha) = @_;
    my $dag_copy = $dag->create_copy();
    my $s_beta = construct($dag_copy, $s_alpha);

    my $considered = {};
    while( !equal($s_beta, $s_alpha) ) {
	my $y = get_rightmost_node($s_alpha, $considered);

	my $swapped;
	do {
	    $swapped = 0;
	    my $z = get_right_neighbor($s_beta, $y);
	    if( defined $z && !occurs_right($s_alpha, $z, $y) ) {
		swap($s_beta, $y, $z);
		if( $dag_copy->has_edge($y, $z) ) {
		    $dag_copy->cover_and_reverse($y, $z);
		}
		$swapped = 1;
	    }
	} while($swapped);
    }
    return $dag_copy;
}


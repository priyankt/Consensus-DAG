#!/usr/bin/perl

use strict;
use warnings;
use CDAG;
use AI::Genetic::Pro;
use Getopt::Long;
use Algorithm::Permute;
use Chart::Gnuplot;

use constant MAX_POPULATION => 500;
use constant INITIAL_SCORE => 10000;

# Input number of dags to be used to get consensus dag
my $num_dags = 3;

# Input number of generations to evolve
my $num_gen = 50;
my $file_prefix = 'dag';
my $all = 0;
my $chart_name = 'chart.gif';
my $random = 1200;
my $sample = 100;

GetOptions(
    "gen=i"    => \$num_gen,
    "chart:s"  => \$chart_name,
    "random:i" => \$random,
    "sample:i" => \$sample,
    "all"      => \$all
) or die "Invalid command line arguments";

my $count = 0;
my $sample_data = {};

my @vertices = (0..9);
my $dags = [];

# Run $sample number of times
while($count < $sample) {

    print STDERR "Processing sample $count\n";
    my $cnt=0;
    $dags = [];
    # Generate 3 random graphs
    while($cnt < 3) {
	my $g = Graph->random_graph(vertices => \@vertices, edges => 12);
	if( $g->is_dag && $g->isolated_vertices() <= 0) {
	    push @$dags, new CDAG($g);
	    $cnt++;
	}
    }

    my $fittest = [];

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
    $ga->evolve($num_gen);

    push @{ $sample_data->{ga} }, $ga->as_value($ga->getFittest);
    push @{ $sample_data->{rand} }, get_random_data_score();

    for(my $i=0; $i < scalar(@$dags); $i++) {
	my $score = fitness( {}, [ $dags->[$i]{graph}->topological_sort ] );
	push @{ $sample_data->{$i} }, $score;
    }
    $count++;
}

my @datasets = ();
my $ga_data = Chart::Gnuplot::DataSet->new(
    xdata => [1..$sample],
    ydata => $sample_data->{ga},
    title => 'GA Best Score',
    style => "lines",
    );
push @datasets, $ga_data;

my $rand_data = Chart::Gnuplot::DataSet->new(
    xdata => [1..$sample],
    ydata => $sample_data->{rand},
    title => 'Random Best Score',
    style => "linespoints",
    );
push @datasets, $rand_data;

my $g0_data = Chart::Gnuplot::DataSet->new(
    xdata => [1..$sample],
    ydata => $sample_data->{0},
    title => 'G0 Score',
    style => "linespoints",
    );
push @datasets, $g0_data;

my $g1_data = Chart::Gnuplot::DataSet->new(
    xdata => [1..$sample],
    ydata => $sample_data->{1},
    title => 'G1 Score',
    style => "linespoints",
    );
push @datasets, $g1_data;

my $g2_data = Chart::Gnuplot::DataSet->new(
    xdata => [1..$sample],
    ydata => $sample_data->{2},
    title => 'G2 Score',
    style => "linespoints",
    );
push @datasets, $g2_data;

my $chart = Chart::Gnuplot->new(
        output => $chart_name,
        title  => "Score Comaprison",
        xlabel => "Iteration",
        ylabel => "Best Score",
    );
$chart->plot2d(@datasets);


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
    
    return $ga->as_value($ga->getFittest);;
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


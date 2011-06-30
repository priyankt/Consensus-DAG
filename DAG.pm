package DAG;

use strict;
use warnings;
use Graph::Directed;

sub new {

    my ($class, $fname) = @_;

    my $graph = new Graph::Directed();

    my $self = {
	graph => $graph,
	edges_added => 0,
	order => [],
    };
    bless $self, $class;

    if($fname) {
	$self->_populate_graph($fname);
	$self->{order} = [ $self->{graph}->topological_sort(empty_if_cyclic => 1) ];
    }

    return $self;
}

sub _populate_graph
{
    my ($self, $fname) = @_;
    open(FH, $fname) or die "Unable to open file - $fname";
    while (<FH>) {
	my ($v1, $v2) = split(',', $_);
	chomp($v1);
	chomp($v2);
	$self->insert_edge($v1, $v2);
    }
    if(!$self->{graph}->is_dag) {
	print STDERR "Given graph is not DAG. Found cycle = ", $self->{graph}->find_a_cycle, "\n";
	exit(1);
    }
    print STDERR "Vertices: ", scalar($self->{graph}->vertices), "\n";
    print STDERR "Edges: ", scalar($self->{graph}->edges), "\n"; 
}

sub to_string
{
    my $self = shift;
    return $self->{graph};
}

sub get_rightmost_sink_node {

    my ($self, $ordering) = @_;
    my $graph = $self->{graph};

    my $sink_node = undef;
    if(@$ordering > 0) {
	foreach my $node (reverse @$ordering) {
	    if( $graph->has_vertex($node) && $graph->is_successorless_vertex($node) ) {
		$sink_node = $node;
		last;
	    }
	}
    }
    return $sink_node;
}

=head2 is_parent()

    Returns true if other_node is parent of current_node

    Params:
    $current_node - Node whose parent is to be checked
    $other_node - Node which is to be check for parent

    Returns:
    true - if $other_node is parent of $current_node
    false - otherwise

=cut
sub is_parent {

    my ($self, $current_node, $other_node) = @_;
    my $graph = $self->{graph};

    return $graph->has_edge($other_node, $current_node);
}

=head2 remove_node()

=cut
sub remove_node {

    my ($self, $node) = @_;
    $self->{graph}->delete_vertex($node);
}

=head2 reverse_edge()
=cut
sub reverse_edge {
    my ($self, $v1, $v2) = @_;
    $self->insert_edge($v2, $v1);
    $self->{graph}->delete_edge($v1, $v2);
}

=head2 create_copy()
    Creates deep copy of current DAG object

=cut
sub create_copy {

    my $self = shift;
    my $copy = $self->{graph}->deep_copy;
    my $self_copy = {
	graph => $copy,
	order => $self->{order},
	edges_added => $self->{edges_added},
    };
    my $class = ref($self);
    bless $self_copy, $class;

    return $self_copy;
}

=head2 insert_node()
=cut
sub insert_node {
    my ($self, $node) = @_;

    my $graph = $self->{graph};
    if( !$graph->has_vertex($node) ) {
	$graph->add_vertex($node);
    }
}

=head2 insert_edge()
=cut
sub insert_edge {
    my ($self, $v1, $v2) = @_;

    my $graph = $self->{graph};
    if( !$graph->has_vertex($v1) ) {
	$graph->add_vertex($v1);
    }
    if(!$graph->has_vertex($v2)) {
	$graph->add_vertex($v2);
    }
    if(!$graph->has_edge($v1, $v2) && $v1 != $v2) {
	$graph->add_edge($v1, $v2);
    }
}

=head2 has_edge()
=cut
sub has_edge {
    my ($self, $v1, $v2) = @_;
    return $self->{graph}->has_edge($v1, $v2);
}

=head2 cover_and_reverse(v1, v2)

    Covers the edge v1->v2 and reverses the edge to v2->v1
    More on 'covering' can be found in paper 'Finding Consensus DAGS'

    Params:
      v1 - parent vertex in edge v1->v2
      v2 - child vertex in edge v1->v2

    Returns:
      None

=cut
sub cover_and_reverse {
    my ($self, $v1, $v2) = @_;

    my $graph = $self->{graph};
    my @parents_v1 = $graph->predecessors($v1);
    my @parents_v2 = $graph->predecessors($v2);
    my %count = ();

    foreach my $item (@parents_v1, @parents_v2 ) {
	$count{$item}++;
    }

    # Cover v1->v2 by adding edges from parents that are missing
    foreach my $item (keys %count) {
	if( $count{$item} == 1 && $item != $v1 ) {
	    $self->insert_edge($item, $v1);
	    $self->insert_edge($item, $v2);
	    $self->{edges_added}++;
	}
    }

    # Reverse v1->v2 to v2->v1
    $self->reverse_edge($v1, $v2);
}

=head2 get_vertices()

    Returns all the vertices in the graph

    Params:
      None

    Returns:
      All vertices in graph - Array/List

=cut
sub get_vertices {
    my $self = shift;
    return $self->{graph}->vertices;
}

=head2 get_cost()

    Returns the total number of edges added to convert graph G to G_alpha

    Params:
      None

    Returns:
      Number of edges added - Int

=cut
sub get_cost {
    my $self = shift;
    return $self->{edges_added};
}

1;

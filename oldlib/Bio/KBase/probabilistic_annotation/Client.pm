package Bio::KBase::probabilistic_annotation::Client;

use JSON::RPC::Client;
use strict;
use Data::Dumper;
use URI;
use Bio::KBase::Exceptions;
use Bio::KBase::AuthToken;

# Client version should match Impl version
# This is a Semantic Version number,
# http://semver.org
our $VERSION = "0.1.0";

=head1 NAME

Bio::KBase::probabilistic_annotation::Client

=head1 DESCRIPTION


The purpose of the Probabilistic Annotation service is to provide users with
alternative annotations for genes, each attached to a likelihood score, and to
translate these likelihood scores into likelihood scores for the existence of
reactions in metabolic models.  With the Probabilistic Annotation service:

- Users can quickly assess the quality of an annotation.

- Reaction likelihood computations allow users to estimate the quality of
  metabolic networks generated using the automated reconstruction tools in
  other services.

- Combining reaction likelihoods with gapfilling both directly incorporates
  available genetic evidence into the gapfilling process and provides putative
  gene annotations automatically, reducing the effort needed to search for
  evidence for gapfilled reactions.


=cut

sub new
{
    my($class, $url, @args) = @_;
    

    my $self = {
	client => Bio::KBase::probabilistic_annotation::Client::RpcClient->new,
	url => $url,
    };

    #
    # This module requires authentication.
    #
    # We create an auth token, passing through the arguments that we were (hopefully) given.

    {
	my $token = Bio::KBase::AuthToken->new(@args);
	
	if (!$token->error_message)
	{
	    $self->{token} = $token->token;
	    $self->{client}->{token} = $token->token;
	}
        else
        {
	    #
	    # All methods in this module require authentication. In this case, if we
	    # don't have a token, we can't continue.
	    #
	    die "Authentication failed: " . $token->error_message;
	}
    }

    my $ua = $self->{client}->ua;	 
    my $timeout = $ENV{CDMI_TIMEOUT} || (30 * 60);	 
    $ua->timeout($timeout);
    bless $self, $class;
    #    $self->_validate_version();
    return $self;
}




=head2 annotate

  $jobid = $obj->annotate($input)

=over 4

=item Parameter and return types

=begin html

<pre>
$input is an AnnotateParams
$jobid is a job_id
AnnotateParams is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	genome_workspace has a value which is a workspace_id
	probanno has a value which is a probanno_id
	probanno_workspace has a value which is a workspace_id
	overwrite has a value which is a bool
	verbose has a value which is a bool
genome_id is a string
workspace_id is a string
probanno_id is a string
bool is an int
job_id is a string

</pre>

=end html

=begin text

$input is an AnnotateParams
$jobid is a job_id
AnnotateParams is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	genome_workspace has a value which is a workspace_id
	probanno has a value which is a probanno_id
	probanno_workspace has a value which is a workspace_id
	overwrite has a value which is a bool
	verbose has a value which is a bool
genome_id is a string
workspace_id is a string
probanno_id is a string
bool is an int
job_id is a string


=end text

=item Description

Generate alternative annotations for every gene in a genome together with
their likelihoods.  Results are stored in a ProbAnno object. Returns the
job ID of the submitted job.

=back

=cut

sub annotate
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function annotate (received $n, expecting 1)");
    }
    {
	my($input) = @args;

	my @_bad_arguments;
        (ref($input) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"input\" (value was \"$input\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to annotate:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'annotate');
	}
    }

    my $result = $self->{client}->call($self->{url}, {
	method => "ProbabilisticAnnotation.annotate",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'annotate',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method annotate",
					    status_line => $self->{client}->status_line,
					    method_name => 'annotate',
				       );
    }
}



=head2 calculate

  $output = $obj->calculate($input)

=over 4

=item Parameter and return types

=begin html

<pre>
$input is a CalculateParams
$output is an object_metadata
CalculateParams is a reference to a hash where the following keys are defined:
	probanno has a value which is a probanno_id
	probanno_workspace has a value which is a workspace_id
	template_model has a value which is a template_id
	template_model_workspace has a value which is a workspace_id
	rxnprobs has a value which is a rxnprobs_id
	rxnprobs_workspace has a value which is a workspace_id
	verbose has a value which is a bool
probanno_id is a string
workspace_id is a string
template_id is a string
rxnprobs_id is a string
bool is an int
object_metadata is a reference to a list containing 11 items:
	0: (id) an object_id
	1: (type) an object_type
	2: (moddate) a timestamp
	3: (instance) an int
	4: (command) a string
	5: (lastmodifier) a username
	6: (owner) a username
	7: (workspace) a workspace_id
	8: (ref) a workspace_ref
	9: (chsum) a string
	10: (metadata) a reference to a hash where the key is a string and the value is a string
object_id is a string
object_type is a string
timestamp is a string
username is a string
workspace_ref is a string

</pre>

=end html

=begin text

$input is a CalculateParams
$output is an object_metadata
CalculateParams is a reference to a hash where the following keys are defined:
	probanno has a value which is a probanno_id
	probanno_workspace has a value which is a workspace_id
	template_model has a value which is a template_id
	template_model_workspace has a value which is a workspace_id
	rxnprobs has a value which is a rxnprobs_id
	rxnprobs_workspace has a value which is a workspace_id
	verbose has a value which is a bool
probanno_id is a string
workspace_id is a string
template_id is a string
rxnprobs_id is a string
bool is an int
object_metadata is a reference to a list containing 11 items:
	0: (id) an object_id
	1: (type) an object_type
	2: (moddate) a timestamp
	3: (instance) an int
	4: (command) a string
	5: (lastmodifier) a username
	6: (owner) a username
	7: (workspace) a workspace_id
	8: (ref) a workspace_ref
	9: (chsum) a string
	10: (metadata) a reference to a hash where the key is a string and the value is a string
object_id is a string
object_type is a string
timestamp is a string
username is a string
workspace_ref is a string


=end text

=item Description

Calculate reaction likelihoods from a probabilistic annotation and a
template model.  Results are stored in a RxnProbs object.  Returns the
metadata for the reaction probability object.

=back

=cut

sub calculate
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function calculate (received $n, expecting 1)");
    }
    {
	my($input) = @args;

	my @_bad_arguments;
        (ref($input) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"input\" (value was \"$input\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to calculate:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'calculate');
	}
    }

    my $result = $self->{client}->call($self->{url}, {
	method => "ProbabilisticAnnotation.calculate",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'calculate',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method calculate",
					    status_line => $self->{client}->status_line,
					    method_name => 'calculate',
				       );
    }
}



=head2 get_rxnprobs

  $output = $obj->get_rxnprobs($input)

=over 4

=item Parameter and return types

=begin html

<pre>
$input is a GetRxnprobsParams
$output is a reaction_probability_list
GetRxnprobsParams is a reference to a hash where the following keys are defined:
	rxnprobs has a value which is a rxnprobs_id
	rxnprobs_workspace has a value which is a workspace_id
rxnprobs_id is a string
workspace_id is a string
reaction_probability_list is a reference to a list where each element is a reaction_probability
reaction_probability is a reference to a list containing 5 items:
	0: (reaction) a reaction_id
	1: (probability) a float
	2: (type) a string
	3: (complex_info) a string
	4: (gene_list) a string
reaction_id is a string

</pre>

=end html

=begin text

$input is a GetRxnprobsParams
$output is a reaction_probability_list
GetRxnprobsParams is a reference to a hash where the following keys are defined:
	rxnprobs has a value which is a rxnprobs_id
	rxnprobs_workspace has a value which is a workspace_id
rxnprobs_id is a string
workspace_id is a string
reaction_probability_list is a reference to a list where each element is a reaction_probability
reaction_probability is a reference to a list containing 5 items:
	0: (reaction) a reaction_id
	1: (probability) a float
	2: (type) a string
	3: (complex_info) a string
	4: (gene_list) a string
reaction_id is a string


=end text

=item Description

Convert a reaction probability object into a human-readable table.

=back

=cut

sub get_rxnprobs
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function get_rxnprobs (received $n, expecting 1)");
    }
    {
	my($input) = @args;

	my @_bad_arguments;
        (ref($input) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"input\" (value was \"$input\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to get_rxnprobs:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'get_rxnprobs');
	}
    }

    my $result = $self->{client}->call($self->{url}, {
	method => "ProbabilisticAnnotation.get_rxnprobs",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'get_rxnprobs',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method get_rxnprobs",
					    status_line => $self->{client}->status_line,
					    method_name => 'get_rxnprobs',
				       );
    }
}



=head2 get_probanno

  $output = $obj->get_probanno($input)

=over 4

=item Parameter and return types

=begin html

<pre>
$input is a GetProbannoParams
$output is a roleset_probabilities
GetProbannoParams is a reference to a hash where the following keys are defined:
	probanno has a value which is a probanno_id
	probanno_workspace has a value which is a workspace_id
probanno_id is a string
workspace_id is a string
roleset_probabilities is a reference to a hash where the key is a feature_id and the value is a reference to a list where each element is a function_probability
feature_id is a string
function_probability is a reference to a list containing 2 items:
	0: (annotation) a string
	1: (probability) a float

</pre>

=end html

=begin text

$input is a GetProbannoParams
$output is a roleset_probabilities
GetProbannoParams is a reference to a hash where the following keys are defined:
	probanno has a value which is a probanno_id
	probanno_workspace has a value which is a workspace_id
probanno_id is a string
workspace_id is a string
roleset_probabilities is a reference to a hash where the key is a feature_id and the value is a reference to a list where each element is a function_probability
feature_id is a string
function_probability is a reference to a list containing 2 items:
	0: (annotation) a string
	1: (probability) a float


=end text

=item Description

Convert a ProbAnno object into a human-readbale table.

=back

=cut

sub get_probanno
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function get_probanno (received $n, expecting 1)");
    }
    {
	my($input) = @args;

	my @_bad_arguments;
        (ref($input) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"input\" (value was \"$input\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to get_probanno:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'get_probanno');
	}
    }

    my $result = $self->{client}->call($self->{url}, {
	method => "ProbabilisticAnnotation.get_probanno",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'get_probanno',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method get_probanno",
					    status_line => $self->{client}->status_line,
					    method_name => 'get_probanno',
				       );
    }
}



sub version {
    my ($self) = @_;
    my $result = $self->{client}->call($self->{url}, {
        method => "ProbabilisticAnnotation.version",
        params => [],
    });
    if ($result) {
        if ($result->is_error) {
            Bio::KBase::Exceptions::JSONRPC->throw(
                error => $result->error_message,
                code => $result->content->{code},
                method_name => 'get_probanno',
            );
        } else {
            return wantarray ? @{$result->result} : $result->result->[0];
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(
            error => "Error invoking method get_probanno",
            status_line => $self->{client}->status_line,
            method_name => 'get_probanno',
        );
    }
}

sub _validate_version {
    my ($self) = @_;
    my $svr_version = $self->version();
    my $client_version = $VERSION;
    my ($cMajor, $cMinor) = split(/\./, $client_version);
    my ($sMajor, $sMinor) = split(/\./, $svr_version);
    if ($sMajor != $cMajor) {
        Bio::KBase::Exceptions::ClientServerIncompatible->throw(
            error => "Major version numbers differ.",
            server_version => $svr_version,
            client_version => $client_version
        );
    }
    if ($sMinor < $cMinor) {
        Bio::KBase::Exceptions::ClientServerIncompatible->throw(
            error => "Client minor version greater than Server minor version.",
            server_version => $svr_version,
            client_version => $client_version
        );
    }
    if ($sMinor > $cMinor) {
        warn "New client version available for Bio::KBase::probabilistic_annotation::Client\n";
    }
    if ($sMajor == 0) {
        warn "Bio::KBase::probabilistic_annotation::Client version is $svr_version. API subject to change.\n";
    }
}

=head1 TYPES



=head2 bool

=over 4



=item Description

*************************************************************************************


=item Definition

=begin html

<pre>
an int
</pre>

=end html

=begin text

an int

=end text

=back



=head2 probanno_id

=over 4



=item Description

A string identifier for a probabilistic annotation object.


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 job_id

=over 4



=item Description

A string identifier for a job object.


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 template_id

=over 4



=item Description

A string identifier for a template model object.


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 rxnprobs_id

=over 4



=item Description

A string identifier for a reaction probabilities object.


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 genome_id

=over 4



=item Description

A string identifier for a genome.


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 feature_id

=over 4



=item Description

A string identifier for a feature.


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 workspace_id

=over 4



=item Description

A string identifier for a workspace. Any string consisting of alphanumeric characters and "-" is acceptable.


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 object_type

=over 4



=item Description

A string indicating the type of an object stored in a workspace.


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 reaction_id

=over 4



=item Description

A string identifier for a reaction object.


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 object_id

=over 4



=item Description

ID of an object stored in the workspace.


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 username

=over 4



=item Description

Login name of KBase user account to which permissions for workspaces are mapped


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 timestamp

=over 4



=item Description

Exact time for workspace operations. e.g. 2012-12-17T23:24:06


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 workspace_ref

=over 4



=item Description

A permanent reference to an object in a workspace.


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 object_metadata

=over 4



=item Description

Meta data associated with an object stored in a workspace.

        object_id id - ID of the object assigned by the user or retreived from the IDserver (e.g. kb|g.0)
        object_type type - type of the object (e.g. Genome)
        timestamp moddate - date when the object was modified by the user (e.g. 2012-12-17T23:24:06)
        int instance - instance of the object, which is equal to the number of times the user has overwritten the object
        timestamp date_created - time at which the alignment was built/loaded in seconds since the epoch
        string command - name of the command last used to modify or create the object
        username lastmodifier - name of the user who last modified the object
        username owner - name of the user who owns (who created) this object
        workspace_id workspace - ID of the workspace in which the object is currently stored
        workspace_ref ref - a 36 character ID that provides permanent undeniable access to this specific instance of this object
        string chsum - checksum of the associated data object
        mapping<string,string> metadata - custom metadata entered for data object during save operation


=item Definition

=begin html

<pre>
a reference to a list containing 11 items:
0: (id) an object_id
1: (type) an object_type
2: (moddate) a timestamp
3: (instance) an int
4: (command) a string
5: (lastmodifier) a username
6: (owner) a username
7: (workspace) a workspace_id
8: (ref) a workspace_ref
9: (chsum) a string
10: (metadata) a reference to a hash where the key is a string and the value is a string

</pre>

=end html

=begin text

a reference to a list containing 11 items:
0: (id) an object_id
1: (type) an object_type
2: (moddate) a timestamp
3: (instance) an int
4: (command) a string
5: (lastmodifier) a username
6: (owner) a username
7: (workspace) a workspace_id
8: (ref) a workspace_ref
9: (chsum) a string
10: (metadata) a reference to a hash where the key is a string and the value is a string


=end text

=back



=head2 function_probability

=over 4



=item Description

A function_probability is a (annotation, probability) pair associated with a gene
An annotation is a "///"-delimited list of roles that could be associated with that gene.


=item Definition

=begin html

<pre>
a reference to a list containing 2 items:
0: (annotation) a string
1: (probability) a float

</pre>

=end html

=begin text

a reference to a list containing 2 items:
0: (annotation) a string
1: (probability) a float


=end text

=back



=head2 ProbAnno

=over 4



=item Description

Object to carry alternative functions and probabilities for genes in a genome    

        probanno_id id - ID of the probabilistic annotation object    
        genome_id genome - ID of the genome the probabilistic annotation was built for
        workspace_id genome_workspace - ID of the workspace containing genome
        mapping<feature_id, list<function_probability>> roleset_probabilities - mapping of features to list of alternative function_probability objects
        list<feature_id> skipped_features - list of features in genome with no probability


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
id has a value which is a probanno_id
genome has a value which is a genome_id
genome_workspace has a value which is a workspace_id
roleset_probabilities has a value which is a reference to a hash where the key is a feature_id and the value is a reference to a list where each element is a function_probability
skipped_features has a value which is a reference to a list where each element is a feature_id

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
id has a value which is a probanno_id
genome has a value which is a genome_id
genome_workspace has a value which is a workspace_id
roleset_probabilities has a value which is a reference to a hash where the key is a feature_id and the value is a reference to a list where each element is a function_probability
skipped_features has a value which is a reference to a list where each element is a feature_id


=end text

=back



=head2 reaction_probability

=over 4



=item Description

Data structure to hold probability of a reaction

        reaction_id reaction - ID of the reaction
        float probability - Probability of the reaction
        string type - Type of complexes ("HASCOMPLEXES" or "NOCOMPLEXES")
        string complex_info - Detailed information on complexes
        string gene_list - List of genes most likely to be attached to reaction


=item Definition

=begin html

<pre>
a reference to a list containing 5 items:
0: (reaction) a reaction_id
1: (probability) a float
2: (type) a string
3: (complex_info) a string
4: (gene_list) a string

</pre>

=end html

=begin text

a reference to a list containing 5 items:
0: (reaction) a reaction_id
1: (probability) a float
2: (type) a string
3: (complex_info) a string
4: (gene_list) a string


=end text

=back



=head2 RxnProbs

=over 4



=item Description

Object to hold reaction probabilities for a genome.

        genome_id genome - ID of the genome the reaction probabilities was built for
        list<reaction_probability> reaction_probabilities - list of reaction probabilities


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
id has a value which is a rxnprobs_id
template_model has a value which is a template_id
template_workspace has a value which is a workspace_id
genome has a value which is a genome_id
genome_workspace has a value which is a workspace_id
probanno has a value which is a probanno_id
probanno_workspace has a value which is a workspace_id
reaction_probabilities has a value which is a reference to a list where each element is a reaction_probability

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
id has a value which is a rxnprobs_id
template_model has a value which is a template_id
template_workspace has a value which is a workspace_id
genome has a value which is a genome_id
genome_workspace has a value which is a workspace_id
probanno has a value which is a probanno_id
probanno_workspace has a value which is a workspace_id
reaction_probabilities has a value which is a reference to a list where each element is a reaction_probability


=end text

=back



=head2 AnnotateParams

=over 4



=item Description

Input parameters for the "annotate" function.

       genome_id genome - ID of Genome object
       workspace_id genome_workspace - ID of workspace where Genome object is stored
       probanno_id probanno - ID of ProbAnno object
       workspace_id probanno_workspace - ID workspace where ProbAnno object is saved
       bool overwrite - True to overwrite existing ProbAnno object with same name
           bool verbose - True to print verbose messages


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
genome has a value which is a genome_id
genome_workspace has a value which is a workspace_id
probanno has a value which is a probanno_id
probanno_workspace has a value which is a workspace_id
overwrite has a value which is a bool
verbose has a value which is a bool

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
genome has a value which is a genome_id
genome_workspace has a value which is a workspace_id
probanno has a value which is a probanno_id
probanno_workspace has a value which is a workspace_id
overwrite has a value which is a bool
verbose has a value which is a bool


=end text

=back



=head2 CalculateParams

=over 4



=item Description

Input parameters for the "calculate" function.

            probanno_id probanno - ID of ProbAnno object
            workspace_id probanno_workspace - ID of workspace where ProbAnno object is stored
            bool verbose - True to print verbose messages


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
probanno has a value which is a probanno_id
probanno_workspace has a value which is a workspace_id
template_model has a value which is a template_id
template_model_workspace has a value which is a workspace_id
rxnprobs has a value which is a rxnprobs_id
rxnprobs_workspace has a value which is a workspace_id
verbose has a value which is a bool

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
probanno has a value which is a probanno_id
probanno_workspace has a value which is a workspace_id
template_model has a value which is a template_id
template_model_workspace has a value which is a workspace_id
rxnprobs has a value which is a rxnprobs_id
rxnprobs_workspace has a value which is a workspace_id
verbose has a value which is a bool


=end text

=back



=head2 GetRxnprobsParams

=over 4



=item Description

Inputs for get_rxnprobs function.
rxnprobs_id rxnprobs- ID for RxnProbs object in the workspace
workspace_id rxnprobs_workspace - ID for workspace in which RxnProbs object is held.


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
rxnprobs has a value which is a rxnprobs_id
rxnprobs_workspace has a value which is a workspace_id

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
rxnprobs has a value which is a rxnprobs_id
rxnprobs_workspace has a value which is a workspace_id


=end text

=back



=head2 reaction_probability_list

=over 4



=item Description

Output for get_rxnprobs function.
It is a list of tuples convenient for output as a table.


=item Definition

=begin html

<pre>
a reference to a list where each element is a reaction_probability
</pre>

=end html

=begin text

a reference to a list where each element is a reaction_probability

=end text

=back



=head2 GetProbannoParams

=over 4



=item Description

Inputs for get_probanno function.
probanno_id probanno - ID for probanno object
workspace_id probanno_workspace - ID for workspace in which ProbAnno object is held.


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
probanno has a value which is a probanno_id
probanno_workspace has a value which is a workspace_id

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
probanno has a value which is a probanno_id
probanno_workspace has a value which is a workspace_id


=end text

=back



=head2 roleset_probabilities

=over 4



=item Description

Output for get_probanno function.
It is a mapping from a feature (gene) ID to a list of (annotation, likelihood) tuples.
Annotations are roles separated by a "///" delimiter


=item Definition

=begin html

<pre>
a reference to a hash where the key is a feature_id and the value is a reference to a list where each element is a function_probability
</pre>

=end html

=begin text

a reference to a hash where the key is a feature_id and the value is a reference to a list where each element is a function_probability

=end text

=back



=cut

package Bio::KBase::probabilistic_annotation::Client::RpcClient;
use base 'JSON::RPC::Client';

#
# Override JSON::RPC::Client::call because it doesn't handle error returns properly.
#

sub call {
    my ($self, $uri, $obj) = @_;
    my $result;

    if ($uri =~ /\?/) {
       $result = $self->_get($uri);
    }
    else {
        Carp::croak "not hashref." unless (ref $obj eq 'HASH');
        $result = $self->_post($uri, $obj);
    }

    my $service = $obj->{method} =~ /^system\./ if ( $obj );

    $self->status_line($result->status_line);

    if ($result->is_success) {

        return unless($result->content); # notification?

        if ($service) {
            return JSON::RPC::ServiceObject->new($result, $self->json);
        }

        return JSON::RPC::ReturnObject->new($result, $self->json);
    }
    elsif ($result->content_type eq 'application/json')
    {
        return JSON::RPC::ReturnObject->new($result, $self->json);
    }
    else {
        return;
    }
}


sub _post {
    my ($self, $uri, $obj) = @_;
    my $json = $self->json;

    $obj->{version} ||= $self->{version} || '1.1';

    if ($obj->{version} eq '1.0') {
        delete $obj->{version};
        if (exists $obj->{id}) {
            $self->id($obj->{id}) if ($obj->{id}); # if undef, it is notification.
        }
        else {
            $obj->{id} = $self->id || ($self->id('JSON::RPC::Client'));
        }
    }
    else {
        # $obj->{id} = $self->id if (defined $self->id);
	# Assign a random number to the id if one hasn't been set
	$obj->{id} = (defined $self->id) ? $self->id : substr(rand(),2);
    }

    my $content = $json->encode($obj);

    $self->ua->post(
        $uri,
        Content_Type   => $self->{content_type},
        Content        => $content,
        Accept         => 'application/json',
	($self->{token} ? (Authorization => $self->{token}) : ()),
    );
}



1;

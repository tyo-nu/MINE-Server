package mineDatabaseServicesImpl;
use strict;
use Bio::KBase::Exceptions;
# Use Semantic Versioning (2.0.0-rc.1)
# http://semver.org 
our $VERSION = "0.1.0";

=head1 NAME

mineDatabaseServices

=head1 DESCRIPTION

=head1 mineDatabaseServices

=head2 SYNOPSIS

The MINE database is fundamentally composed of two different types of documents, which are represented by the Compound
and Reaction objects. Users can use text-matching queries to access these records directly or perform two types of more
advanced queries: Mass Adduct queries and pathway queries. Mass Adduct queries return a list of compounds that might
match the m/z of an unknown compound. Pathway queries return either the shortest path or all paths between two compounds
 in the database.

=cut

#BEGIN_HEADER
#END_HEADER

sub new
{
    my($class, @args) = @_;
    my $self = {
    };
    bless $self, $class;
    #BEGIN_CONSTRUCTOR
    #END_CONSTRUCTOR

    if ($self->can('_init_instance'))
    {
	$self->_init_instance();
    }
    return $self;
}

=head1 METHODS



=head2 quick_search

  $quick_search_results = $obj->quick_search($db, $query)

=over 4

=item Parameter and return types

=begin html

<pre>
$db is a string
$query is a string
$quick_search_results is a reference to a list where each element is a comp_stub
comp_stub is a reference to a hash where the following keys are defined:
	id has a value which is a string
	Formula has a value which is a string
	Mass has a value which is a float
	Inchi_key has a value which is a string

</pre>

=end html

=begin text

$db is a string
$query is a string
$quick_search_results is a reference to a list where each element is a comp_stub
comp_stub is a reference to a hash where the following keys are defined:
	id has a value which is a string
	Formula has a value which is a string
	Mass has a value which is a float
	Inchi_key has a value which is a string


=end text



=item Description

Creates quick_search_results, a list of comp_stubs which match the query string. Searches for matches to KEGG
Codes, Inchi Keys, Brenda IDs and Names.

=back

=cut

sub quick_search
{
    my $self = shift;
    my($db, $query) = @_;

    my @_bad_arguments;
    (!ref($db)) or push(@_bad_arguments, "Invalid type for argument \"db\" (value was \"$db\")");
    (!ref($query)) or push(@_bad_arguments, "Invalid type for argument \"query\" (value was \"$query\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to quick_search:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'quick_search');
    }

    my $ctx = $mineDatabaseServicesServer::CallContext;
    my($quick_search_results);
    #BEGIN quick_search
    #END quick_search
    my @_bad_returns;
    (ref($quick_search_results) eq 'ARRAY') or push(@_bad_returns, "Invalid type for return variable \"quick_search_results\" (value was \"$quick_search_results\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to quick_search:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'quick_search');
    }
    return($quick_search_results);
}




=head2 quick_search

  $quick_search_results = $obj->quick_search($db, $query)

=over 4

=item Parameter and return types

=begin html

<pre>
$db is a string
$query is a string
$quick_search_results is a reference to a list where each element is a comp_stub
comp_stub is a reference to a hash where the following keys are defined:
	id has a value which is a string
	Formula has a value which is a string
	Mass has a value which is a float
	Inchi_key has a value which is a string

</pre>

=end html

=begin text

$db is a string
$query is a string
$quick_search_results is a reference to a list where each element is a comp_stub
comp_stub is a reference to a hash where the following keys are defined:
	id has a value which is a string
	Formula has a value which is a string
	Mass has a value which is a float
	Inchi_key has a value which is a string


=end text



=item Description

Creates similarity_search_results, a list of comp_stubs whose Tannimoto coefficient to the search smiles is
greater that the user set threshold. Uses open babel FP2 fingerprints to match.

=back

=cut

sub quick_search
{
    my $self = shift;
    my($db, $query) = @_;

    my @_bad_arguments;
    (!ref($db)) or push(@_bad_arguments, "Invalid type for argument \"db\" (value was \"$db\")");
    (!ref($query)) or push(@_bad_arguments, "Invalid type for argument \"query\" (value was \"$query\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to quick_search:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'quick_search');
    }

    my $ctx = $mineDatabaseServicesServer::CallContext;
    my($quick_search_results);
    #BEGIN quick_search
    #END quick_search
    my @_bad_returns;
    (ref($quick_search_results) eq 'ARRAY') or push(@_bad_returns, "Invalid type for return variable \"quick_search_results\" (value was \"$quick_search_results\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to quick_search:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'quick_search');
    }
    return($quick_search_results);
}




=head2 database_query

  $database_query_results = $obj->database_query($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a database_query_params
$database_query_results is a reference to a list where each element is a comp_stub
database_query_params is a reference to a hash where the following keys are defined:
	db has a value which is a string
	field has a value which is a string
	value has a value which is a string
	regex has a value which is a bool
bool is an int
comp_stub is a reference to a hash where the following keys are defined:
	id has a value which is a string
	Formula has a value which is a string
	Mass has a value which is a float
	Inchi_key has a value which is a string

</pre>

=end html

=begin text

$params is a database_query_params
$database_query_results is a reference to a list where each element is a comp_stub
database_query_params is a reference to a hash where the following keys are defined:
	db has a value which is a string
	field has a value which is a string
	value has a value which is a string
	regex has a value which is a bool
bool is an int
comp_stub is a reference to a hash where the following keys are defined:
	id has a value which is a string
	Formula has a value which is a string
	Mass has a value which is a float
	Inchi_key has a value which is a string


=end text



=item Description

Creates database_query_results, a list of object_ids which match the json query string

=back

=cut

sub database_query
{
    my $self = shift;
    my($params) = @_;

    my @_bad_arguments;
    (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"params\" (value was \"$params\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to database_query:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'database_query');
    }

    my $ctx = $mineDatabaseServicesServer::CallContext;
    my($database_query_results);
    #BEGIN database_query
    #END database_query
    my @_bad_returns;
    (ref($database_query_results) eq 'ARRAY') or push(@_bad_returns, "Invalid type for return variable \"database_query_results\" (value was \"$database_query_results\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to database_query:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'database_query');
    }
    return($database_query_results);
}




=head2 get_comps

  $objects = $obj->get_comps($db, $ids)

=over 4

=item Parameter and return types

=begin html

<pre>
$db is a string
$ids is a reference to a list where each element is an object_id
$objects is a reference to a list where each element is a CompoundObject
object_id is a string
CompoundObject is a reference to a hash where the following keys are defined:
	id has a value which is an object_id
	InChI_Key has a value which is a string
	Formula has a value which is a string
	Stringcode has a value which is a string
	Mass has a value which is a float
	Charge has a value which is an int
	KEGG_Code has a value which is a reference to a list where each element is a string
	BRENDA_Name has a value which is a reference to a list where each element is a string
	Reactant_in has a value which is a reference to a list where each element is an object_id
	Product_of has a value which is a reference to a list where each element is an object_id

</pre>

=end html

=begin text

$db is a string
$ids is a reference to a list where each element is an object_id
$objects is a reference to a list where each element is a CompoundObject
object_id is a string
CompoundObject is a reference to a hash where the following keys are defined:
	id has a value which is an object_id
	InChI_Key has a value which is a string
	Formula has a value which is a string
	Stringcode has a value which is a string
	Mass has a value which is a float
	Charge has a value which is an int
	KEGG_Code has a value which is a reference to a list where each element is a string
	BRENDA_Name has a value which is a reference to a list where each element is a string
	Reactant_in has a value which is a reference to a list where each element is an object_id
	Product_of has a value which is a reference to a list where each element is an object_id


=end text



=item Description

Return a list of CompoundObjects that match supplied object_ids in a specified db

=back

=cut

sub get_comps
{
    my $self = shift;
    my($db, $ids) = @_;

    my @_bad_arguments;
    (!ref($db)) or push(@_bad_arguments, "Invalid type for argument \"db\" (value was \"$db\")");
    (ref($ids) eq 'ARRAY') or push(@_bad_arguments, "Invalid type for argument \"ids\" (value was \"$ids\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to get_comps:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'get_comps');
    }

    my $ctx = $mineDatabaseServicesServer::CallContext;
    my($objects);
    #BEGIN get_comps
    #END get_comps
    my @_bad_returns;
    (ref($objects) eq 'ARRAY') or push(@_bad_returns, "Invalid type for return variable \"objects\" (value was \"$objects\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to get_comps:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'get_comps');
    }
    return($objects);
}




=head2 get_rxns

  $objects = $obj->get_rxns($db, $ids)

=over 4

=item Parameter and return types

=begin html

<pre>
$db is a string
$ids is a reference to a list where each element is an object_id
$objects is a reference to a list where each element is a ReactionObject
object_id is a string
ReactionObject is a reference to a hash where the following keys are defined:
	id has a value which is an object_id
	Operators has a value which is a reference to a list where each element is a string
	Reactants has a value which is a reference to a list where each element is a rxn_comp
	Products has a value which is a reference to a list where each element is a rxn_comp
	Energy has a value which is a float
	Error has a value which is a float
rxn_comp is a reference to a list containing 2 items:
	0: (stoic) an int
	1: (id) an object_id

</pre>

=end html

=begin text

$db is a string
$ids is a reference to a list where each element is an object_id
$objects is a reference to a list where each element is a ReactionObject
object_id is a string
ReactionObject is a reference to a hash where the following keys are defined:
	id has a value which is an object_id
	Operators has a value which is a reference to a list where each element is a string
	Reactants has a value which is a reference to a list where each element is a rxn_comp
	Products has a value which is a reference to a list where each element is a rxn_comp
	Energy has a value which is a float
	Error has a value which is a float
rxn_comp is a reference to a list containing 2 items:
	0: (stoic) an int
	1: (id) an object_id


=end text



=item Description

Returns a list of ReactionObjects that match supplied object_ids in a specified db

=back

=cut

sub get_rxns
{
    my $self = shift;
    my($db, $ids) = @_;

    my @_bad_arguments;
    (!ref($db)) or push(@_bad_arguments, "Invalid type for argument \"db\" (value was \"$db\")");
    (ref($ids) eq 'ARRAY') or push(@_bad_arguments, "Invalid type for argument \"ids\" (value was \"$ids\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to get_rxns:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'get_rxns');
    }

    my $ctx = $mineDatabaseServicesServer::CallContext;
    my($objects);
    #BEGIN get_rxns
    #END get_rxns
    my @_bad_returns;
    (ref($objects) eq 'ARRAY') or push(@_bad_returns, "Invalid type for return variable \"objects\" (value was \"$objects\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to get_rxns:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'get_rxns');
    }
    return($objects);
}




=head2 get_models

  $models = $obj->get_models()

=over 4

=item Parameter and return types

=begin html

<pre>
$models is a reference to a list where each element is a reference to a list containing 2 items:
	0: (id) a string
	1: (name) a string

</pre>

=end html

=begin text

$models is a reference to a list where each element is a reference to a list containing 2 items:
	0: (id) a string
	1: (name) a string


=end text



=item Description

Returns a list of SEED models available to be set as native metabolites as tuples of SEED id and name

=back

=cut

sub get_models
{
    my $self = shift;

    my $ctx = $mineDatabaseServicesServer::CallContext;
    my($models);
    #BEGIN get_models
    #END get_models
    my @_bad_returns;
    (ref($models) eq 'ARRAY') or push(@_bad_returns, "Invalid type for return variable \"models\" (value was \"$models\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to get_models:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'get_models');
    }
    return($models);
}




=head2 get_adducts

  $adducts = $obj->get_adducts()

=over 4

=item Parameter and return types

=begin html

<pre>
$adducts is a reference to a list containing 2 items:
	0: a reference to a list where each element is a string
	1: a reference to a list where each element is a string

</pre>

=end html

=begin text

$adducts is a reference to a list containing 2 items:
	0: a reference to a list where each element is a string
	1: a reference to a list where each element is a string


=end text



=item Description

Returns a tuple of lists of positive and negative mass adducts names that may be used for querying the databases

=back

=cut

sub get_adducts
{
    my $self = shift;

    my $ctx = $mineDatabaseServicesServer::CallContext;
    my($adducts);
    #BEGIN get_adducts
    #END get_adducts
    my @_bad_returns;
    (ref($adducts) eq 'ARRAY') or push(@_bad_returns, "Invalid type for return variable \"adducts\" (value was \"$adducts\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to get_adducts:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'get_adducts');
    }
    return($adducts);
}




=head2 adduct_db_search

  $output = $obj->adduct_db_search($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a mass_adduct_query_params
$output is a reference to a list where each element is an adduct_result
mass_adduct_query_params is a reference to a hash where the following keys are defined:
	db has a value which is a string
	mz has a value which is a float
	tolerance has a value which is a float
	adduct_list has a value which is a reference to a list where each element is a string
	models has a value which is a reference to a list where each element is a string
	ppm has a value which is a bool
	charge has a value which is a bool
	halogens has a value which is a bool
bool is an int
adduct_result is a reference to a hash where the following keys are defined:
	adduct has a value which is a string
	formula has a value which is a string
	isomers has a value which is a reference to a list where each element is an object_id
object_id is a string

</pre>

=end html

=begin text

$params is a mass_adduct_query_params
$output is a reference to a list where each element is an adduct_result
mass_adduct_query_params is a reference to a hash where the following keys are defined:
	db has a value which is a string
	mz has a value which is a float
	tolerance has a value which is a float
	adduct_list has a value which is a reference to a list where each element is a string
	models has a value which is a reference to a list where each element is a string
	ppm has a value which is a bool
	charge has a value which is a bool
	halogens has a value which is a bool
bool is an int
adduct_result is a reference to a hash where the following keys are defined:
	adduct has a value which is a string
	formula has a value which is a string
	isomers has a value which is a reference to a list where each element is an object_id
object_id is a string


=end text



=item Description

Creates output, a list of adduct, formula and isomer combinations that match the supplied parameters

=back

=cut

sub adduct_db_search
{
    my $self = shift;
    my($params) = @_;

    my @_bad_arguments;
    (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"params\" (value was \"$params\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to adduct_db_search:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'adduct_db_search');
    }

    my $ctx = $mineDatabaseServicesServer::CallContext;
    my($output);
    #BEGIN adduct_db_search
    #END adduct_db_search
    my @_bad_returns;
    (ref($output) eq 'ARRAY') or push(@_bad_returns, "Invalid type for return variable \"output\" (value was \"$output\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to adduct_db_search:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'adduct_db_search');
    }
    return($output);
}




=head2 pathway_search

  $pathway_query_results = $obj->pathway_search($pathway_query_params)

=over 4

=item Parameter and return types

=begin html

<pre>
$pathway_query_params is a pathway_query_params
$pathway_query_results is a reference to a list where each element is a pathway
pathway_query_params is a reference to a hash where the following keys are defined:
	db has a value which is a string
	start_comp has a value which is an object_id
	end_comp has a value which is an object_id
	len_limit has a value which is an int
	all_paths has a value which is a bool
	np_min has a value which is a float
	gibbs_cap has a value which is a float
object_id is a string
bool is an int
pathway is a reference to a list where each element is an object_id

</pre>

=end html

=begin text

$pathway_query_params is a pathway_query_params
$pathway_query_results is a reference to a list where each element is a pathway
pathway_query_params is a reference to a hash where the following keys are defined:
	db has a value which is a string
	start_comp has a value which is an object_id
	end_comp has a value which is an object_id
	len_limit has a value which is an int
	all_paths has a value which is a bool
	np_min has a value which is a float
	gibbs_cap has a value which is a float
object_id is a string
bool is an int
pathway is a reference to a list where each element is an object_id


=end text



=item Description

Creates pathway_query_results, a list of valid pathways (length one unless all_paths is true)

=back

=cut

sub pathway_search
{
    my $self = shift;
    my($pathway_query_params) = @_;

    my @_bad_arguments;
    (ref($pathway_query_params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"pathway_query_params\" (value was \"$pathway_query_params\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to pathway_search:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'pathway_search');
    }

    my $ctx = $mineDatabaseServicesServer::CallContext;
    my($pathway_query_results);
    #BEGIN pathway_search
    #END pathway_search
    my @_bad_returns;
    (ref($pathway_query_results) eq 'ARRAY') or push(@_bad_returns, "Invalid type for return variable \"pathway_query_results\" (value was \"$pathway_query_results\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to pathway_search:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'pathway_search');
    }
    return($pathway_query_results);
}




=head2 version 

  $return = $obj->version()

=over 4

=item Parameter and return types

=begin html

<pre>
$return is a string
</pre>

=end html

=begin text

$return is a string

=end text

=item Description

Return the module version. This is a Semantic Versioning number.

=back

=cut

sub version {
    return $VERSION;
}

=head1 TYPES



=head2 bool

=over 4



=item Description

indicates true or false values, false = 0, true =1


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



=head2 object_id

=over 4



=item Description

Unique ID of a compound or reaction derived from a hexdigest of the sha1 hash of a unique feature.
Starts with C if a compound, X if a cofactor and R if a reaction.


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



=head2 comp_stub

=over 4



=item Description

A summery of a compound object which is returned from compound query

        string _id - unique ID of a compound
        string Formula - molecular formula of the compound
        float Mass - exact mass of the compound
        string Inchi_key - the Inchi Key of the compound


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
id has a value which is a string
Formula has a value which is a string
Mass has a value which is a float
Inchi_key has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
id has a value which is a string
Formula has a value which is a string
Mass has a value which is a float
Inchi_key has a value which is a string


=end text

=back



=head2 rxn_comp

=over 4



=item Description

A compound that is a component of a reaction as tuple of stoichiometric coefficient and _id


=item Definition

=begin html

<pre>
a reference to a list containing 2 items:
0: (stoic) an int
1: (id) an object_id

</pre>

=end html

=begin text

a reference to a list containing 2 items:
0: (stoic) an int
1: (id) an object_id


=end text

=back



=head2 pathway

=over 4



=item Description

A list of all the compounds and reactions in a pathway


=item Definition

=begin html

<pre>
a reference to a list where each element is an object_id
</pre>

=end html

=begin text

a reference to a list where each element is an object_id

=end text

=back



=head2 peak

=over 4



=item Description

An annotated ms peak output by a batch mass adduct query(not yet implemented)

        string name - name of the peak
        int num_forms - number of formula hits
        int num_hits - total number of compound matches


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
name has a value which is a string
num_forms has a value which is an int
num_hits has a value which is an int

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
name has a value which is a string
num_forms has a value which is an int
num_hits has a value which is an int


=end text

=back



=head2 adduct_result

=over 4



=item Description

The result of a single adduct query on the database

        string adduct - the name of the mass adduct that returned the result
        string formula; - the formula that was matched
        list<object_id> - a list of the isomers of the formula present in the database


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
adduct has a value which is a string
formula has a value which is a string
isomers has a value which is a reference to a list where each element is an object_id

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
adduct has a value which is a string
formula has a value which is a string
isomers has a value which is a reference to a list where each element is an object_id


=end text

=back



=head2 CompoundObject

=over 4



=item Description

Data structures for a compound object

                Guaranteed:
                object_id id - A hexdigest of the sha1 hash of the openbabel canonical smile
                string InChI_Key - The first block of the InChI Key of a compound
                string Formula - The chemical formula of the compound
        string Stringcode - The canonical SMILE string generated by openbabel
                float Mass - The exact mass of the neutral form of a compound as calculated by openbabel
                int Charge - The total charge of the compound as calculated by ChemAxon

                Optionally:
                list<string> KEGG_Code - KEGG compound codes
        list<string> BRENDA_Name - Names from the BRENDA repository
        list<object_id> Reactant_in - Reactions in which the compound is a reactant
        list<object_id> Product_of - Reactions in which the compound is a product


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
id has a value which is an object_id
InChI_Key has a value which is a string
Formula has a value which is a string
Stringcode has a value which is a string
Mass has a value which is a float
Charge has a value which is an int
KEGG_Code has a value which is a reference to a list where each element is a string
BRENDA_Name has a value which is a reference to a list where each element is a string
Reactant_in has a value which is a reference to a list where each element is an object_id
Product_of has a value which is a reference to a list where each element is an object_id

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
id has a value which is an object_id
InChI_Key has a value which is a string
Formula has a value which is a string
Stringcode has a value which is a string
Mass has a value which is a float
Charge has a value which is an int
KEGG_Code has a value which is a reference to a list where each element is a string
BRENDA_Name has a value which is a reference to a list where each element is a string
Reactant_in has a value which is a reference to a list where each element is an object_id
Product_of has a value which is a reference to a list where each element is an object_id


=end text

=back



=head2 ReactionObject

=over 4



=item Description

Data structures for a reaction object

                Guaranteed:
                object_id id - A hexdigest of the sha1 hash of the _ids of the reactants and products in sorted order
        list<string> Operators - The operator used to generate a particular reaction
        rxn_comps Reactants - Reactants of the reaction as tuples
        rxn_comps Products - Products of the reaction as tuples

        Optionally:
        float Energy - Delta G of reaction calculated by group contribution theory
        float Error - Estimated error of above energy


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
id has a value which is an object_id
Operators has a value which is a reference to a list where each element is a string
Reactants has a value which is a reference to a list where each element is a rxn_comp
Products has a value which is a reference to a list where each element is a rxn_comp
Energy has a value which is a float
Error has a value which is a float

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
id has a value which is an object_id
Operators has a value which is a reference to a list where each element is a string
Reactants has a value which is a reference to a list where each element is a rxn_comp
Products has a value which is a reference to a list where each element is a rxn_comp
Energy has a value which is a float
Error has a value which is a float


=end text

=back



=head2 database_query_params

=over 4



=item Description

Input parameters for the "database_query" function.

        string db - the database against which the query will be performed
string field - the field of the database to match
string value - the value to match
bool regex - if true the value will be processed as a regular expression


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
db has a value which is a string
field has a value which is a string
value has a value which is a string
regex has a value which is a bool

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
db has a value which is a string
field has a value which is a string
value has a value which is a string
regex has a value which is a bool


=end text

=back



=head2 mass_adduct_query_params

=over 4



=item Description

Input parameters for the "mass_adduct_query" function.

        string db - the database in which to search for mass spec matches
        float mz - the experimental mass per charge ratio
float tolerance - the desired mass precision
list<adduct> adduct_list - the adducts to consider in the query.
list<string> models - the models in SEED that will be considered native metabolites
string charge - the polarity for molecules if not specified by file
bool ppm - if true, precision is supplied in parts per million. Else, precision is in Daltons
bool halogens - if false, compounds containing Cl, Br, and F will be excluded from results


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
db has a value which is a string
mz has a value which is a float
tolerance has a value which is a float
adduct_list has a value which is a reference to a list where each element is a string
models has a value which is a reference to a list where each element is a string
ppm has a value which is a bool
charge has a value which is a bool
halogens has a value which is a bool

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
db has a value which is a string
mz has a value which is a float
tolerance has a value which is a float
adduct_list has a value which is a reference to a list where each element is a string
models has a value which is a reference to a list where each element is a string
ppm has a value which is a bool
charge has a value which is a bool
halogens has a value which is a bool


=end text

=back



=head2 pathway_query_params

=over 4



=item Description

Input parameters for the "pathway_search" function.

        string db - the database in which to search for pathways
        object_id start_comp - the compound to begin the search from
object_id end_comp - the compound that that a pathway will end with if successful
int len_limit - the max number of intermediate reactions permitted in a path.
bool all_paths - if true, the script returns all paths less that the limit not just the shortest path
float np_min - Set a floor on the minimum natural product likeness of any one compound in a pathway
float gibbs_cap - Set a cap on the gibbs free energy of any one reaction in a pathway


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
db has a value which is a string
start_comp has a value which is an object_id
end_comp has a value which is an object_id
len_limit has a value which is an int
all_paths has a value which is a bool
np_min has a value which is a float
gibbs_cap has a value which is a float

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
db has a value which is a string
start_comp has a value which is an object_id
end_comp has a value which is an object_id
len_limit has a value which is an int
all_paths has a value which is a bool
np_min has a value which is a float
gibbs_cap has a value which is a float


=end text

=back



=cut

1;

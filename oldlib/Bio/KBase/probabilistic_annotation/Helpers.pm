package Bio::KBase::probabilistic_annotation::Helpers;

# Need a license here

use strict;
use warnings;
use Exporter;
use parent qw(Exporter);
our @EXPORT_OK = qw( get_probanno_url get_probanno_client );
our $defaultURL = "http://kbase.us/services/probabilistic_annotation/";

# Get the URL for the Probabilistic Annotation service.
sub get_probanno_url {
    my $set = shift;
    my $currentURL;
    if (defined($set)) {
    	if ($set eq "default") {
            $set = $defaultURL;
        }
    	$currentURL = $set;
		my $filename;
    	if (!defined($ENV{KB_RUNNING_IN_IRIS})) {
			$filename = "$ENV{HOME}/.kbase_probannoURL";
    	} else {
			$filename = "/.kbase_probannoURL";
    	}
	    open(my $fh, ">", $filename) || return;
		print $fh $currentURL;
		close($fh);
    } elsif (!defined($currentURL)) {
		my $filename;
    	if (!defined($ENV{KB_RUNNING_IN_IRIS})) {
		    $filename = "$ENV{HOME}/.kbase_probannoURL";
		} else {
			$filename = "/.kbase_probannoURL";
		}
	    if	( -e $filename ) {
			open(my $fh, "<", $filename) || return;
			$currentURL = <$fh>;
			chomp $currentURL;
			close($fh);
		} else {
			$currentURL = $defaultURL;
		}
    }
    return $currentURL;
}

# Return a new Probabilistic Annotation client using the currently set URL.
sub get_probanno_client() {
	return Bio::KBase::probabilistic_annotation::Client->new(get_probanno_url());
}
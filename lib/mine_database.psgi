use mineDatabaseServicesImpl;

use mineDatabaseServicesServer;
use Plack::Middleware::CrossOrigin;



my @dispatch;

{
    my $obj = mineDatabaseServicesImpl->new;
    push(@dispatch, 'mineDatabaseServices' => $obj);
}


my $server = mineDatabaseServicesServer->new(instance_dispatch => { @dispatch },
				allow_get => 0,
			       );

my $handler = sub { $server->handle_input(@_) };

$handler = Plack::Middleware::CrossOrigin->wrap( $handler, origins => "*", headers => "*");

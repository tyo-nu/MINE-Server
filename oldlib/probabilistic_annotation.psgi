use ProbabilisticAnnotationImpl;

use ProbabilisticAnnotationServer;
use Plack::Middleware::CrossOrigin;



my @dispatch;

{
    my $obj = ProbabilisticAnnotationImpl->new;
    push(@dispatch, 'ProbabilisticAnnotation' => $obj);
}


my $server = ProbabilisticAnnotationServer->new(instance_dispatch => { @dispatch },
				allow_get => 0,
			       );

my $handler = sub { $server->handle_input(@_) };

$handler = Plack::Middleware::CrossOrigin->wrap( $handler, origins => "*", headers => "*");

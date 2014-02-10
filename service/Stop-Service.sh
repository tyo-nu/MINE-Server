pidfile=/vol/model-prod/mine-server/service/pid2
uwsgi --stop $pidfile
rm $pidfile

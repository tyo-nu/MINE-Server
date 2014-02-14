pidfile=/vol/model-prod/mine-server/service/pid
uwsgi --stop $pidfile
rm $pidfile

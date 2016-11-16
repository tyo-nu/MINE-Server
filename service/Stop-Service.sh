source /vol/kbase/deployment/user-env.sh
pidfile=/vol/model-prod/mine-server/service/pid
uwsgi --stop $pidfile
rm $pidfile
